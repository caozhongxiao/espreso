
#include "generatedinput.h"

#include "../basis/containers/serializededata.h"
#include "../basis/logging/logging.h"
#include "../basis/logging/timeeval.h"
#include "../basis/utilities/communication.h"
#include "../basis/utilities/utils.h"
#include "../config/ecf/environment.h"

#include "../mesh/mesh.h"
#include "../mesh/store/nodestore.h"

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
	ESINFO(OVERVIEW) << "Build mesh from generated data.";
	TimeEval timing("Load generated mesh");
	timing.totalTime.startWithBarrier();

	TimeEvent tdangling("remove dangling nodes"); tdangling.start();
	removeDanglingNodes();
	tdangling.end(); timing.addEvent(tdangling);
	ESINFO(PROGRESS2) << "Generated data loader:: dangling nodes removed.";

	TimeEvent tneighs("fill neighbors"); tneighs.start();
	fillNeighbors();
	tneighs.end(); timing.addEvent(tneighs);
	ESINFO(PROGRESS2) << "Generated data loader:: neighbors filled.";

	if (needSynchronization) {
		TimeEvent tsync("synchronize global indices"); tsync.start();
		synchronizeGlobalIndices();
		tsync.end(); timing.addEvent(tsync);
		ESINFO(PROGRESS2) << "Generated data loader:: global indices synchronized.";

		TimeEvent tnsort("sort nodes"); tnsort.start();
		sortNodes(true);
		tnsort.end(); timing.addEvent(tnsort);
		ESINFO(PROGRESS2) << "Generated data loader:: nodes sorted.";
	}

	TimeEvent tnodes("fill nodes"); tnodes.start();
	fillNodes();
	tnodes.end(); timing.addEvent(tnodes);
	ESINFO(PROGRESS2) << "Generated data loader:: nodes filled.";

	TimeEvent telements("fill elements"); telements.start();
	fillElements();
	telements.end(); timing.addEvent(telements);
	ESINFO(PROGRESS2) << "Generated data loader:: elements filled.";

	TimeEvent tregions("fill regions"); tregions.start();
	fillElementRegions();
	fillBoundaryRegions();
	fillNodeRegions();
	tregions.end(); timing.addEvent(tregions);
	ESINFO(PROGRESS2) << "Generated data loader:: elements filled.";

	timing.totalTime.endWithBarrier();
	timing.printStatsMPI();
}

void GeneratedInput::removeDanglingNodes()
{
	std::vector<eslocal> usedNodes = _meshData.enodes;
	Esutils::sortAndRemoveDuplicity(usedNodes);

	std::vector<Point> coordinates;
	std::vector<eslocal> nIDs, ndist, noffset(_meshData.nIDs.size());
	std::vector<int> nranks;

	ndist.push_back(0);
	for (size_t i = 0; i < usedNodes.size(); i++) {
		noffset[usedNodes[i]] = i;
		coordinates.push_back(_meshData.coordinates[usedNodes[i]]);
		nIDs.push_back(_meshData.nIDs[usedNodes[i]]);
		nranks.insert(nranks.end(), _meshData.nranks.begin() + _meshData.ndist[usedNodes[i]], _meshData.nranks.begin() + _meshData.ndist[usedNodes[i] + 1]);
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
	_meshData.ndist.swap(ndist);
	_meshData.nranks.swap(nranks);
}

struct __Point__ {

	static constexpr size_t digits = 1000;

	__Point__(): x(0), y(0), z(0), id(-1) {};
	__Point__(const Point &p, esglobal id): x(p.x), y(p.y), z(p.z), id(id) {};

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
	esglobal id;
};

void GeneratedInput::synchronizeGlobalIndices()
{
	auto n2i = [ & ] (size_t neighbour) {
		return std::lower_bound(_mesh.neighbours.begin(), _mesh.neighbours.end(), neighbour) - _mesh.neighbours.begin();
	};

	std::vector<std::vector<__Point__> > sBuffer(_mesh.neighbours.size());
	std::vector<std::vector<__Point__> > rBuffer(_mesh.neighbours.size());

	for (size_t n = 0; n < _meshData.nIDs.size(); ++n) {
		if (_meshData.ndist[n + 1] - _meshData.ndist[n] > 1) {
			if (_meshData.nranks[_meshData.ndist[n]] == environment->MPIrank) {
				for (eslocal r = _meshData.ndist[n] + 1; r < _meshData.ndist[n + 1]; ++r) {
					sBuffer[n2i(_meshData.nranks[r])].push_back(__Point__(_meshData.coordinates[n], _meshData.nIDs[n]));
				}
			} else {
				sBuffer[n2i(_meshData.nranks[_meshData.ndist[n]])].push_back(__Point__(_meshData.coordinates[n], n));
			}
		}
	}

	for (size_t n = 0; n < sBuffer.size(); n++) {
		std::sort(sBuffer[n].begin(), sBuffer[n].end());
		rBuffer[n].resize(sBuffer[n].size());
	}

	if (!Communication::receiveLowerKnownSize(sBuffer, rBuffer, _mesh.neighbours)) {
		ESINFO(ERROR) << "problem while synchronization of global indices.";
	}

	for (size_t n = 0; n < _mesh.neighbours.size(); n++) {
		if (_mesh.neighbours[n] < environment->MPIrank) {
			for (size_t p = 0; p < rBuffer[n].size(); p++) {
				auto it = std::lower_bound(sBuffer[n].begin(), sBuffer[n].end(), rBuffer[n][p]);
				if (*it == rBuffer[n][p]) {
					_meshData.nIDs[it->id] = rBuffer[n][p].id;
				} else {
					ESINFO(ERROR) << "Internal ERROR while synchronization global indices: " << _mesh.neighbours[n] << " on " << environment->MPIrank;
				}
			}
		}
	}
}


