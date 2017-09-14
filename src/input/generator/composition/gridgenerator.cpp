
#include "gridgenerator.h"

#include <random>

#include "../../../config/ecf/environment.h"
#include "../../../config/ecf/input/grid.h"
#include "../../../basis/utilities/parser.h"
#include "../primitives/block.h"
#include "../generator.h"
#include "../../../mesh/structures/region.h"
#include "../../../mesh/structures/mesh.h"
#include "../../../mesh/structures/elementtypes.h"
#include "../../../mesh/structures/coordinates.h"

#include "../elements/2D/square4.h"
#include "../elements/2D/square8.h"
#include "../elements/2D/triangle3.h"
#include "../elements/2D/triangle6.h"

#include "../elements/3D/hexahedron20.h"
#include "../elements/3D/hexahedron8.h"
#include "../elements/3D/prisma15.h"
#include "../elements/3D/prisma6.h"
#include "../elements/3D/pyramid13.h"
#include "../elements/3D/pyramid5.h"
#include "../elements/3D/tetrahedron10.h"
#include "../elements/3D/tetrahedron4.h"

using namespace espreso::input;

GridSettings::GridSettings()
: etype(GENERATOR_ELEMENT_TYPE::HEXA8),
  start(0, 0, 0), end(1, 1, 1),
  blocks(1, 1, 1), clusters(1, 1, 1), domains(2, 2, 2), elements(5, 5, 5),
  projection(Expression("x", { "x", "y", "z" }), Expression("y", { "x", "y", "z" }), Expression("z", { "x", "y", "z" })),
  rotation(Expression("0", { "x", "y", "z" }), Expression("0", { "x", "y", "z" }), Expression("0", { "x", "y", "z" })),
  nonempty(1, true), uniformDecomposition(true) {}

GridSettings::GridSettings(const GridGeneratorConfiguration &configuration)
: projection(Expression(configuration.projection_x, { "x", "y", "z" }), Expression(configuration.projection_y, { "x", "y", "z" }), Expression(configuration.projection_z, { "x", "y", "z" })),
  rotation(Expression(configuration.rotation_x, { "x", "y", "z" }), Expression(configuration.rotation_y, { "x", "y", "z" }), Expression(configuration.rotation_z, { "x", "y", "z" }))
{
	etype = configuration.element_type;

	blocks     = Triple<size_t>(configuration.blocks_x, configuration.blocks_y, configuration.blocks_z);
	clusters = Triple<size_t>(configuration.clusters_x, configuration.clusters_y, configuration.clusters_z);
	nonempty.resize((blocks * clusters).mul(), true);

	domains  = Triple<size_t>(configuration.domains_x, configuration.domains_y, configuration.domains_z);
	elements = Triple<size_t>(configuration.elements_x, configuration.elements_y, configuration.elements_z);

	start = Triple<esglobal>(configuration.start_x / Generator::precision, configuration.start_y / Generator::precision, configuration.start_z / Generator::precision);
	end   = Triple<esglobal>(
			(configuration.start_x + configuration.length_x) / Generator::precision,
			(configuration.start_y + configuration.length_y) / Generator::precision,
			(configuration.start_z + configuration.length_z) / Generator::precision);

	for (auto it = configuration.blocks.begin(); it != configuration.blocks.end(); ++it) {
		if (it->first >= nonempty.size()) {
			ESINFO(GLOBAL_ERROR) << "Block index is out of range.";
		}
		nonempty[it->first] = it->second;
	}

	uniformDecomposition = configuration.uniform_decomposition;
}

void GridGenerator::load(const GridGeneratorConfiguration &configuration, Mesh &mesh, size_t index, size_t size)
{
	GridGenerator grid(configuration, mesh, index, size);
	grid.fill();
}

GridGenerator::GridGenerator(const GridGeneratorConfiguration &configuration, Mesh &mesh, size_t index, size_t size)
: Loader(mesh), _configuration(configuration), _settings(configuration), _index(index), _size(size), _body(0)
{
	Triple<size_t> clusters = _settings.blocks * _settings.clusters;
	std::string element;

	size_t cluster = 0;
	Triple<size_t> offset;
	int clusterIndex = 0;
	for (offset.z = 0; offset.z < clusters.z; offset.z++) {
		for (offset.y = 0; offset.y < clusters.y; offset.y++) {
			for (offset.x = 0; offset.x < clusters.x; offset.x++, clusterIndex++) {

				if (_settings.nonempty[cluster]) {
					_cMap.push_back(cluster);
				} else {
					_cMap.push_back(-1);
					continue;
				}

				if (cluster++ == index) {
					_clusterOffset = offset;
					_clusterIndex = clusterIndex;
					BlockSetting block;
					block.domains  = _settings.domains;
					block.elements = _settings.elements;
					block.start = _settings.start + ((_settings.end - _settings.start) / (Triple<double>)_settings.clusters * offset).round();
					block.end   = _settings.start + ((_settings.end - _settings.start) / (Triple<double>)_settings.clusters * (offset + 1)).round();
					block.projection = _settings.projection;
					block.rotation = _settings.rotation;
					switch (_settings.etype) {
					case GENERATOR_ELEMENT_TYPE::HEXA8:
						_block = new Block<Hexahedron8>(block);
						_subnodes = Hexahedron8::subnodes;
						element = "hexahedrons";
						break;
					case GENERATOR_ELEMENT_TYPE::HEXA20:
						_block = new Block<Hexahedron20>(block);
						_subnodes = Hexahedron20::subnodes;
						element = "hexahedrons with midnodes";
						break;
					case GENERATOR_ELEMENT_TYPE::TETRA4:
						_block = new Block<Tetrahedron4>(block);
						_subnodes = Tetrahedron4::subnodes;
						element = "tetrahedrons";
						break;
					case GENERATOR_ELEMENT_TYPE::TETRA10:
						_block = new Block<Tetrahedron10>(block);
						_subnodes = Tetrahedron10::subnodes;
						element = "tetrahedrons with midnodes";
						break;
					case GENERATOR_ELEMENT_TYPE::PRISMA6:
						_block = new Block<Prisma6>(block);
						_subnodes = Prisma6::subnodes;
						element = "prismas";
						break;
					case GENERATOR_ELEMENT_TYPE::PRISMA15:
						_block = new Block<Prisma15>(block);
						_subnodes = Prisma15::subnodes;
						element = "prismas with midnodes";
						break;
					case GENERATOR_ELEMENT_TYPE::PYRAMID5:
						_block = new Block<Pyramid5>(block);
						_subnodes = Pyramid5::subnodes;
						element = "pyramids";
						break;
					case GENERATOR_ELEMENT_TYPE::PYRAMID13:
						_block = new Block<Pyramid13>(block);
						_subnodes = Pyramid13::subnodes;
						element = "pyramids with midnodes";
						break;

					case GENERATOR_ELEMENT_TYPE::SQUARE4:
						_block = new Block<Square4>(block);
						_subnodes = Square4::subnodes;
						element = "squares";
						break;
					case GENERATOR_ELEMENT_TYPE::SQUARE8:
						_block = new Block<Square8>(block);
						_subnodes = Square8::subnodes;
						element = "squares with midnodes";
						break;
					case GENERATOR_ELEMENT_TYPE::TRIANGLE3:
						_block = new Block<Triangle3>(block);
						_subnodes = Triangle3::subnodes;
						element = "triangles";
						break;
					case GENERATOR_ELEMENT_TYPE::TRIANGLE6:
						_block = new Block<Triangle6>(block);
						_subnodes = Triangle6::subnodes;
						element = "triangles with midnodes";
						break;
					}
				}
			}
		}
	}
	if (cluster != size) {
		ESINFO(GLOBAL_ERROR) << "Incorrect number of MPI processes (" << size << "). Should be " << cluster;
	}
	ESINFO(OVERVIEW) << "Generate grid of " << element;
}

GridGenerator::~GridGenerator()
{
	delete _block;
}

size_t GridGenerator::pointCount() const
{
	return (_settings.blocks * _settings.clusters * _settings.domains * _settings.elements * (Triple<size_t>(_subnodes) - 1) + 1).mul();
}

void GridGenerator::points(Coordinates &coordinates, size_t globalIdOffset)
{
	_block->points(coordinates._points);

	Triple<size_t> cnodes = _settings.domains * _settings.elements * (Triple<size_t>(_subnodes) - 1);
	Triple<size_t> gnodes = _settings.blocks * _settings.clusters * cnodes;
	Triple<size_t> coffset = _clusterOffset * cnodes;
	Triple<size_t> size = (gnodes + 1).toSize();

	Triple<size_t> offset;
	for (offset.z = 0; offset.z <= cnodes.z; offset.z++) {
		for (offset.y = 0; offset.y <= cnodes.y; offset.y++) {
			for (offset.x = 0; offset.x <= cnodes.x; offset.x++) {
				coordinates._globalIndex.push_back(((coffset + offset) * size).sum() + globalIdOffset);
				coordinates._globalMapping.push_back(G2L(coordinates._globalIndex.back(), (esglobal)coordinates._globalMapping.size()));
			}
		}
	}
}

void GridGenerator::elements(std::vector<size_t> &bodies, std::vector<Element*> &elements, std::vector<Element*> &faces, std::vector<Element*> &edges)
{
	_block->elements(elements, _body);
	bodies = { 0, elements.size() };
}

bool GridGenerator::partitiate(const std::vector<Element*> &nodes, std::vector<eslocal> &partsPtrs, std::vector<std::vector<Element*> > &fixPoints, std::vector<Element*> &corners)
{
	size_t parts = _settings.domains.mul();
	if (_settings.uniformDecomposition) {
		_block->uniformPartition(partsPtrs, parts);
		_block->uniformFixPoints(nodes, fixPoints);
		_block->uniformCorners(nodes, corners, 1, true, true, true);
		return true;
	} else {
		if (_configuration.random_partition) {
			std::random_device rd;
			std::mt19937 gen(rd());
			std::uniform_int_distribution<> dis(0, parts);
			parts += dis(gen) - parts / 2;
		}
		if (_configuration.noncontinuous.find(_clusterIndex) != _configuration.noncontinuous.end()) {
			mesh.partitiateNoncontinuously(parts, _configuration.noncontinuous.find(_clusterIndex)->second);
		} else {
			mesh.partitiate(parts);
		}
		return false;
	}
}

void GridGenerator::neighbours(std::vector<Element*> &nodes, std::vector<int> &neighbours, const std::vector<Element*> &faces, const std::vector<Element*> &edges)
{
	std::vector<int> map(27);

	Triple<int> offset, coffset(_clusterOffset), count = _settings.blocks * _settings.clusters, size = count.toSize();
	size_t index = 0;
	for (offset.z = -1; offset.z <= 1; offset.z++) {
		for (offset.y = -1; offset.y <= 1; offset.y++) {
			for (offset.x = -1; offset.x <= 1; offset.x++, index++) {
				if ((coffset + offset) < 0 || (count - (coffset + offset)) < 1) {
					map[index] = -1;
				} else {
					map[index] = _cMap[((coffset + offset) * size).sum()];
					if (map[index] != environment->MPIrank) {
						neighbours.push_back(map[index]);
					}
				}
			}
		}
	}

	_block->boundaries(nodes, map);
	std::sort(neighbours.begin(), neighbours.end());
}

void GridGenerator::regions(
		std::vector<Evaluator*> &evaluators,
		std::vector<Region*> &regions,
		std::vector<Element*> &elements,
		std::vector<Element*> &faces,
		std::vector<Element*> &edges,
		std::vector<Element*> &nodes)
{
	for (auto it = _configuration.nodes.begin(); it != _configuration.nodes.end(); ++it) {
		if (StringCompare::caseInsensitiveEq("all", it->second)) {
			regions.push_back(new Region(ElementType::NODES, nodes));
		} else {
			regions.push_back(new Region(ElementType::NODES));
			BlockBorder border(it->second);
			_block->region(nodes, regions.back(), border, 0);
		}
		regions.back()->name = it->first;
	}
	for (auto it = _configuration.edges.begin(); it != _configuration.edges.end(); ++it) {
		if (StringCompare::caseInsensitiveEq("all", it->second)) {
			regions.push_back(new Region(ElementType::EDGES));
			ESINFO(GLOBAL_ERROR) << "Implement region of all edges.";
		} else {
			regions.push_back(new Region(ElementType::EDGES));
			BlockBorder border(it->second);
			_block->region(nodes, regions.back(), border, 1);
			edges.insert(edges.end(), regions.back()->elements().begin(), regions.back()->elements().end());
		}
		regions.back()->name = it->first;
	}
	for (auto it = _configuration.faces.begin(); it != _configuration.faces.end(); ++it) {
		if (StringCompare::caseInsensitiveEq("all", it->second)) {
			regions.push_back(new Region(ElementType::FACES));
			ESINFO(GLOBAL_ERROR) << "Implement region of all faces.";
		} else {
			regions.push_back(new Region(ElementType::FACES));
			BlockBorder border(it->second);
			_block->region(nodes, regions.back(), border, 2);
			faces.insert(faces.end(), regions.back()->elements().begin(), regions.back()->elements().end());
		}
		regions.back()->name = it->first;
	}

	for (auto it = _configuration.elements.begin(); it != _configuration.elements.end(); ++it) {
		if (StringCompare::caseInsensitiveEq("all", it->second)) {
			regions.push_back(new Region(ElementType::ELEMENTS, elements));
		} else if (StringCompare::caseInsensitiveEq("chessboard_white", it->second)) {
			regions.push_back(new Region(ElementType::ELEMENTS));
			_block->pattern(elements, regions.back(), _clusterOffset, _settings.blocks * _settings.clusters, Pattern::CHESSBOARD_WHITE, _configuration.chessboard_size);
		} else if (StringCompare::caseInsensitiveEq("chessboard_black", it->second)) {
			regions.push_back(new Region(ElementType::ELEMENTS));
			_block->pattern(elements, regions.back(), _clusterOffset, _settings.blocks * _settings.clusters, Pattern::CHESSBOARD_BLACK, _configuration.chessboard_size);
		} else if (StringCompare::caseInsensitiveEq("not_selected", it->second)) {
			// skip and process after all regions are set
			continue;
		} else {
			regions.push_back(new Region(ElementType::ELEMENTS));
			BlockBorder border(it->second);
			_block->region(elements, regions.back(), border, 3);
		}
		regions.back()->name = it->first;
	}
	for (auto it = _configuration.elements.begin(); it != _configuration.elements.end(); ++it) {
		if (StringCompare::caseInsensitiveEq("not_selected", it->second)) {
			std::vector<Element*> all(mesh.elements()), selected;
			for (size_t r = 2; r < regions.size(); r++) {
				if (regions[r]->eType == ElementType::ELEMENTS) {
					selected.insert(selected.end(), regions[r]->elements().begin(), regions[r]->elements().end());
				}
			}
			std::sort(all.begin(), all.end());
			std::sort(selected.begin(), selected.end());
			Esutils::removeDuplicity(selected);

			if (selected.size()) {
				regions.push_back(new Region(ElementType::ELEMENTS));
				for (size_t e = 0, sindex = 0; e < all.size(); e++) {
					if (all[e] == selected[sindex]) {
						sindex++;
					} else {
						regions.back()->elements().push_back(all[e]);
					}
				}
			} else {
				regions.push_back(new Region(ElementType::ELEMENTS, elements));
			}
			regions.back()->name = it->first;
		}
	}
}