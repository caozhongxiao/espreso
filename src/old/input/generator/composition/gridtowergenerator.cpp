
#include "gridtowergenerator.h"

#include "../../../config/ecf/environment.h"
#include "../../../config/ecf/input/gridtower.h"

#include "../../../old/mesh/structures/mesh.h"
#include "../../../old/mesh/structures/elementtypes.h"
#include "../primitives/block.h"
#include "../generator.h"

#include "mpi.h"

using namespace espreso::input;

void GridTowerGenerator::load(const GridTowerGeneratorConfiguration &configuration, Mesh &mesh, size_t index, size_t size)
{
	GridTowerGenerator gridtower(configuration, mesh, index, size);
	gridtower.fill();
}

GridTowerGenerator::GridTowerGenerator(const GridTowerGeneratorConfiguration &configuration, Mesh &mesh, size_t index, size_t size)
: OldLoader(mesh), _configuration(configuration), _gridGenerator(NULL), _gridPointsIDOffset(0), _index(index), _size(size)
{
	size_t gridIndex = 0, firstCluster = 0, lastCluster = 0;
	Triple<esglobal> gridPointsOffset(0, 0, 0);
	for (auto it = _configuration.grids.begin(); it != _configuration.grids.end(); ++it, gridIndex++) {
		lastCluster += it->second.clusters_x * it->second.clusters_y * it->second.clusters_z;
		if (firstCluster <= _index && _index < lastCluster) {
			_gridGenerator = new GridGenerator(it->second, mesh, _index - firstCluster, lastCluster - firstCluster);
			_clusterIndexBegin = firstCluster;
			_clusterIndexEnd = lastCluster;
			_gridIndex = gridIndex;
			_gridPointsOffset = gridPointsOffset;
		}
		firstCluster = lastCluster;
		gridPointsOffset += Triple<esglobal>(it->second.length_x / Generator::precision, it->second.length_y / Generator::precision, it->second.length_z / Generator::precision);
	}

	if (lastCluster != size) {
		ESINFO(GLOBAL_ERROR) << "Incorrect number of MPI processes (" << size << "). Should be " << lastCluster;
	}

	size_t gridPointsCounter = _clusterIndexBegin == index ? _gridGenerator->pointCount() : 0;

	std::vector<size_t> counters(size);
	MPI_Allgather(&gridPointsCounter, sizeof(size_t), MPI_BYTE, counters.data(), sizeof(size_t), MPI_BYTE, environment->MPICommunicator);
	for (size_t i = 0; i < _clusterIndexBegin; i++) {
		_gridPointsIDOffset += counters[i];
	}
	_gridGenerator->bodyIndex(_gridIndex);
}

GridTowerGenerator::~GridTowerGenerator()
{
	delete _gridGenerator;
}

void GridTowerGenerator::points(Coordinates &coordinates)
{
	Triple<double> end = _gridGenerator->_block->block.end, start = _gridGenerator->_block->block.start;
	switch (_configuration.direction) {
	case GridTowerGeneratorConfiguration::DIRECTION::X:
		_gridGenerator->_block->block.start.x += _gridPointsOffset.x;
		_gridGenerator->_block->block.end.x += _gridPointsOffset.x;
		break;
	case GridTowerGeneratorConfiguration::DIRECTION::Y:
		_gridGenerator->_block->block.start.y += _gridPointsOffset.y;
		_gridGenerator->_block->block.end.y += _gridPointsOffset.y;
		break;
	case GridTowerGeneratorConfiguration::DIRECTION::Z:
		_gridGenerator->_block->block.start.z += _gridPointsOffset.z;
		_gridGenerator->_block->block.end.z += _gridPointsOffset.z;
		break;
	}

	_gridGenerator->points(coordinates, _gridPointsIDOffset);

	_gridGenerator->_block->block.end = end;
	_gridGenerator->_block->block.start = start;
}

void GridTowerGenerator::elements(std::vector<size_t> &bodies, std::vector<Element*> &elements, std::vector<Element*> &faces, std::vector<Element*> &edges)
{
	_gridGenerator->elements(bodies, elements, faces, edges);
	bodies.resize(_configuration.grids.size() + 1);
	for (size_t i = 0; i < bodies.size(); i++) {
		if (i <= _gridIndex) {
			bodies[i] = 0;
		} else {
			bodies[i] = elements.size();
		}
	}
}

bool GridTowerGenerator::partitiate(const std::vector<Element*> &nodes, std::vector<eslocal> &partsPtrs, std::vector<std::vector<Element*> > &fixPoints, std::vector<Element*> &corners)
{
	return _gridGenerator->partitiate(nodes, partsPtrs, fixPoints, corners);
}

void GridTowerGenerator::neighbours(std::vector<Element*> &nodes, std::vector<int> &neighbours, const std::vector<Element*> &faces, const std::vector<Element*> &edges)
{
	auto middle = _configuration.grids.find(_gridIndex);
	auto lower  = _configuration.grids.end();
	auto upper  = _configuration.grids.end();

	switch (_configuration.direction) {
	case GridTowerGeneratorConfiguration::DIRECTION::X:
		if (_gridGenerator->_clusterOffset.x == 0) {
			lower  = _configuration.grids.find(_gridIndex - 1);
		}
		if (_gridGenerator->_clusterOffset.x == _gridGenerator->_settings.clusters.x - 1) {
			upper  = _configuration.grids.find(_gridIndex + 1);
		}
		break;
	case GridTowerGeneratorConfiguration::DIRECTION::Y:
		if (_gridGenerator->_clusterOffset.y == 0) {
			lower  = _configuration.grids.find(_gridIndex - 1);
		}
		if (_gridGenerator->_clusterOffset.y == _gridGenerator->_settings.clusters.y - 1) {
			upper  = _configuration.grids.find(_gridIndex + 1);
		}
		break;
	case GridTowerGeneratorConfiguration::DIRECTION::Z:
		if (_gridGenerator->_clusterOffset.z == 0) {
			lower  = _configuration.grids.find(_gridIndex - 1);
		}
		if (_gridGenerator->_clusterOffset.z == _gridGenerator->_settings.clusters.z - 1) {
			upper  = _configuration.grids.find(_gridIndex + 1);
		}
		break;
	}

	if (lower != _configuration.grids.end()) {
		if (lower->second.blocks_z != middle->second.blocks_z || lower->second.blocks_y != middle->second.blocks_y || lower->second.blocks_x != middle->second.blocks_x) {
			lower = _configuration.grids.end();
		}
		if (lower->second.clusters_z != middle->second.clusters_z || lower->second.clusters_y != middle->second.clusters_y || lower->second.clusters_x != middle->second.clusters_x) {
			lower = _configuration.grids.end();
		}
		if (lower->second.domains_z != middle->second.domains_z || lower->second.domains_y != middle->second.domains_y || lower->second.domains_x != middle->second.domains_x) {
			lower = _configuration.grids.end();
		}
		if (lower->second.elements_z != middle->second.elements_z || lower->second.elements_y != middle->second.elements_y || lower->second.elements_x != middle->second.elements_x) {
			lower = _configuration.grids.end();
		}
	}
	if (upper != _configuration.grids.end()) {
		if (upper->second.blocks_z != middle->second.blocks_z || upper->second.blocks_y != middle->second.blocks_y || upper->second.blocks_x != middle->second.blocks_x) {
			upper = _configuration.grids.end();
		}
		if (upper->second.clusters_z != middle->second.clusters_z || upper->second.clusters_y != middle->second.clusters_y || upper->second.clusters_x != middle->second.clusters_x) {
			upper = _configuration.grids.end();
		}
		if (upper->second.domains_z != middle->second.domains_z || upper->second.domains_y != middle->second.domains_y || upper->second.domains_x != middle->second.domains_x) {
			upper = _configuration.grids.end();
		}
		if (upper->second.elements_z != middle->second.elements_z || upper->second.elements_y != middle->second.elements_y || upper->second.elements_x != middle->second.elements_x) {
			upper = _configuration.grids.end();
		}
	}

	std::vector<int> map(27);

	Triple<int> offset, coffset(_gridGenerator->_clusterOffset), count = _gridGenerator->_settings.blocks * _gridGenerator->_settings.clusters, size = count.toSize();
	size_t index = 0;
	for (offset.z = -1; offset.z <= 1; offset.z++) {
		for (offset.y = -1; offset.y <= 1; offset.y++) {
			for (offset.x = -1; offset.x <= 1; offset.x++, index++) {
				if ((coffset + offset) < 0 || (count - (coffset + offset)) < 1) {
					map[index] = -1;
				} else {
					map[index] = _gridGenerator->_cMap[((coffset + offset) * size).sum()] + _clusterIndexBegin;
					if (map[index] != environment->MPIrank) {
						neighbours.push_back(map[index]);
					}
				}
			}
		}
	}


	size_t cx = middle->second.blocks_x * middle->second.clusters_x;
	size_t cy = middle->second.blocks_y * middle->second.clusters_y;
	size_t cz = middle->second.blocks_z * middle->second.clusters_z;

	for (offset.z = -1, index = 0; offset.z <= 1; offset.z++) {
		for (offset.y = -1; offset.y <= 1; offset.y++) {
			for (offset.x = -1; offset.x <= 1; offset.x++, index++) {
				switch (_configuration.direction) {
				case GridTowerGeneratorConfiguration::DIRECTION::X:
					if (offset.x == 1 && map[index - 1] != -1 && upper != _configuration.grids.end()) {
						map[index] = map[index - 1] + cx * cy * (cz - 1) + cx * (cy - 1) + 1;
						neighbours.push_back(map[index]);
					}
					if (offset.x == -1 && map[index + 1] != -1 && lower != _configuration.grids.end()) {
						map[index] = map[index + 1] - cx * cy * (cz - 1) - cx * (cy - 1) - 1;
						neighbours.push_back(map[index]);
					}
					break;
				case GridTowerGeneratorConfiguration::DIRECTION::Y:
					if (offset.y == 1 && map[index - 3] != -1 && upper != _configuration.grids.end()) {
						map[index] = map[index - 3] + cx * cy * (cz - 1) + cx;
						neighbours.push_back(map[index]);
					}
					if (offset.y == -1 && map[index + 3] != -1 && lower != _configuration.grids.end()) {
						map[index] = map[index + 3] - cx * cy * (cz - 1) - cx;
						neighbours.push_back(map[index]);
					}
					break;
				case GridTowerGeneratorConfiguration::DIRECTION::Z:
					if (offset.z == 1 && map[index - 9] != -1 && upper != _configuration.grids.end()) {
						map[index] = map[index - 9] + cx * cy;
						neighbours.push_back(map[index]);
					}
					if (offset.z == -1 && map[index + 9] != -1 && lower != _configuration.grids.end()) {
						map[index] = map[index + 9] - cx * cy;
						neighbours.push_back(map[index]);
					}
					break;
				}
			}
		}
	}

	_gridGenerator->_block->boundaries(nodes, map);
	std::sort(neighbours.begin(), neighbours.end());
	mesh.synchronizeGlobalIndices();
}

void GridTowerGenerator::regions(
		std::vector<Evaluator*> &evaluators,
		std::vector<Region*> &regions,
		std::vector<Element*> &elements,
		std::vector<Element*> &faces,
		std::vector<Element*> &edges,
		std::vector<Element*> &nodes)
{
	_gridGenerator->regions(evaluators, regions, elements, faces, edges, nodes);

	auto addRegions = [&] (const std::map<std::string, std::string> &interval, ElementType type) {
		for (auto r = interval.begin(); r != interval.end(); ++r) {
			auto it = std::find_if(regions.begin(), regions.end(), [&] (const Region *region) { return region->name.compare(r->first) == 0; });
			if (it == mesh.regions().end()) {
				regions.push_back(new Region(type));
				regions.back()->name = r->first;
			}
		}
	};


	for (auto it = _configuration.grids.begin(); it != _configuration.grids.end(); ++it) {
		addRegions(it->second.nodes, ElementType::NODES);
		addRegions(it->second.edges, ElementType::EDGES);
		addRegions(it->second.faces, ElementType::FACES);
		addRegions(it->second.elements, ElementType::ELEMENTS);
	}
	mesh.synchronizeRegionOrder();
}



