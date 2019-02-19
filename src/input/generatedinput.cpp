
#include "generatedinput.h"

#include "basis/containers/serializededata.h"
#include "basis/utilities/communication.h"
#include "basis/utilities/utils.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.hpp"

#include "mesh/mesh.h"
#include "mesh/store/nodestore.h"

#include <algorithm>
#include <numeric>


using namespace espreso;

void GeneratedInput::buildMesh(PlainMeshData &meshData, Mesh &mesh, bool needSynchronization)
{
	GeneratedInput(meshData, mesh, needSynchronization);
}

GeneratedInput::GeneratedInput(PlainMeshData &meshData, Mesh &mesh, bool needSynchronization)
: Input(meshData, mesh)
{
	eslog::startln("MESIO: BUILD GENERATED MESH");

	removeDanglingNodes();
	eslog::checkpointln("MESIO: DANGLING NODES REMOVED");

	fillNeighbors();
	eslog::checkpointln("MESIO: NEIGHBOURS FILLED");

	if (needSynchronization) {
		synchronizeGlobalIndices();
		eslog::checkpointln("MESIO: NODES INDICES SYNCHRONIZED");

		sortNodes(true);
		eslog::checkpointln("MESIO: NODES SORTED");
	}

	fillNodes();
	eslog::checkpointln("MESIO: NODES FILLED");

	fillElements();
	eslog::checkpointln("MESIO: ELEMENTS FILLED");

	fillElementRegions();
	fillBoundaryRegions();
	fillNodeRegions();
	eslog::endln("MESIO: REGIONS FILLED");
}

void GeneratedInput::removeDanglingNodes()
{
	std::vector<esint> usedNodes = _meshData.enodes;
	utils::sortAndRemoveDuplicity(usedNodes);

	std::vector<Point> coordinates;
	std::vector<esint> nIDs, ndist, noffset(_meshData.nIDs.size());
	std::vector<int> nranks;

	ndist.push_back(0);
	for (size_t i = 0; i < usedNodes.size(); i++) {
		noffset[usedNodes[i]] = i;
		coordinates.push_back(_meshData.coordinates[usedNodes[i]]);
		nIDs.push_back(_meshData.nIDs[usedNodes[i]]);
		nranks.insert(nranks.end(), _meshData._nranks.begin() + _meshData._nrankdist[usedNodes[i]], _meshData._nranks.begin() + _meshData._nrankdist[usedNodes[i] + 1]);
		ndist.push_back(nranks.size());
	}

	for (size_t i = 0; i < _meshData.enodes.size(); i++) {
		_meshData.enodes[i] = noffset[_meshData.enodes[i]];
	}

	for (auto it = _meshData.nregions.begin(); it != _meshData.nregions.end(); ++it) {
		for (size_t i = 0; i < it->second.size(); i++) {
			it->second[i] = noffset[it->second[i]];
		}
	}

	_meshData.nIDs.swap(nIDs);
	_meshData.coordinates.swap(coordinates);
	_meshData._nrankdist.swap(ndist);
	_meshData._nranks.swap(nranks);
}

struct __Point__ {

	static constexpr size_t digits = 1000;

	__Point__(): x(0), y(0), z(0), id(-1) {};
	__Point__(const Point &p, esint id): x(p.x), y(p.y), z(p.z), id(id) {};

	bool operator<(const __Point__ &p) const
	{
		if (std::trunc(z * digits) == std::trunc(p.z * digits)) {
			if (std::trunc(y * digits) == std::trunc(p.y * digits)) {
				if (std::trunc(x * digits) == std::trunc(p.x * digits)) {
					return false;
				}
				return x < p.x;
			}
			return y < p.y;
		}
		return z < p.z;
	}

	bool operator==(const __Point__ &p) const
	{
		return !(*this < p) && !(p < *this);
	}

	double x, y, z;
	esint id;
};

void GeneratedInput::synchronizeGlobalIndices()
{
	auto n2i = [ & ] (size_t neighbour) {
		return std::lower_bound(_mesh.neighbours.begin(), _mesh.neighbours.end(), neighbour) - _mesh.neighbours.begin();
	};

	std::vector<std::vector<__Point__> > sBuffer(_mesh.neighbours.size());
	std::vector<std::vector<__Point__> > rBuffer(_mesh.neighbours.size());

	for (size_t n = 0; n < _meshData.nIDs.size(); ++n) {
		if (_meshData._nrankdist[n + 1] - _meshData._nrankdist[n] > 1) {
			if (_meshData._nranks[_meshData._nrankdist[n]] == info::mpi::rank) {
				for (esint r = _meshData._nrankdist[n] + 1; r < _meshData._nrankdist[n + 1]; ++r) {
					sBuffer[n2i(_meshData._nranks[r])].push_back(__Point__(_meshData.coordinates[n], _meshData.nIDs[n]));
				}
			} else {
				sBuffer[n2i(_meshData._nranks[_meshData._nrankdist[n]])].push_back(__Point__(_meshData.coordinates[n], n));
			}
		}
	}

	for (size_t n = 0; n < sBuffer.size(); n++) {
		std::sort(sBuffer[n].begin(), sBuffer[n].end());
		rBuffer[n].resize(sBuffer[n].size());
	}

	if (!Communication::receiveLowerKnownSize(sBuffer, rBuffer, _mesh.neighbours)) {
		eslog::error("problem while synchronization of global indices.\n");
	}

	for (size_t n = 0; n < _mesh.neighbours.size(); n++) {
		if (_mesh.neighbours[n] < info::mpi::rank) {
			for (size_t p = 0; p < rBuffer[n].size(); p++) {
				auto it = std::lower_bound(sBuffer[n].begin(), sBuffer[n].end(), rBuffer[n][p]);
				if (*it == rBuffer[n][p]) {
					_meshData.nIDs[it->id] = rBuffer[n][p].id;
				} else {
					eslog::error("Internal ERROR while synchronization global indices.\n");
				}
			}
		}
	}
}


