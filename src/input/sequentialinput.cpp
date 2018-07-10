
#include "sequentialinput.h"

#include "../basis/containers/serializededata.h"
#include "../basis/logging/logging.h"
#include "../basis/logging/timeeval.h"
#include "../basis/utilities/utils.h"

#include "../config/ecf/environment.h"

#include "../mesh/mesh.h"
#include "../mesh/store/nodestore.h"
#include "../mesh/store/elementstore.h"
#include "../mesh/store/boundaryregionstore.h"

#include <algorithm>
#include <numeric>

using namespace espreso;

void SequentialInput::buildMesh(const ECFRoot &configuration, PlainMeshData &meshData, Mesh &mesh)
{
	SequentialInput(configuration, meshData, mesh);
}

SequentialInput::SequentialInput(const ECFRoot &configuration, PlainMeshData &meshData, Mesh &mesh)
: Input(configuration, meshData, mesh)
{
	ESINFO(OVERVIEW) << "Build mesh.";
	TimeEval timing("Load sequential mesh");
	timing.totalTime.startWithBarrier();

	TimeEvent tnsort("sort nodes"); tnsort.start();
	_nDistribution = { 0, _meshData.nIDs.size() };
	sortNodes();
	tnsort.end(); timing.addEvent(tnsort);
	ESINFO(PROGRESS2) << "Sequential loader:: nodes sorted.";

	TimeEvent tesort("sort elements"); tesort.start();
	_eDistribution = { 0, _meshData.eIDs.size() };
	sortElements();
	tesort.end(); timing.addEvent(tesort);
	ESINFO(PROGRESS2) << "Sequential loader:: elements sorted.";

	TimeEvent tranks("fill ranks"); tranks.start();
	_meshData.nranks.resize(_meshData.nIDs.size());
	_meshData.ndist.resize(_meshData.nIDs.size() + 1);
	std::iota(_meshData.ndist.begin(), _meshData.ndist.end(), 0);
	_mesh.neighboursWithMe.push_back(0);
	tranks.end(); timing.addEvent(tranks);
	ESINFO(PROGRESS2) << "Sequential loader:: ranks filled.";

	TimeEvent tpost("regions prepared"); tpost.start();
	prepareRegions();
	tpost.end(); timing.addEvent(tpost);
	ESINFO(PROGRESS2) << "Sequential loader:: regions prepared.";

	TimeEvent tnodes("fill nodes"); tnodes.start();
	fillSortedNodes();
	tnodes.end(); timing.addEvent(tnodes);
	ESINFO(PROGRESS2) << "Sequential loader:: nodes filled.";

	TimeEvent telements("fill elements"); telements.start();
	fillSortedElements();
	telements.end(); timing.addEvent(telements);
	ESINFO(PROGRESS2) << "Sequential loader:: elements filled.";

	TimeEvent treindex("reindex elements nodes"); treindex.start();
	if (
			!std::is_sorted(_mesh.nodes->IDs->datatarray().begin(), _mesh.nodes->IDs->datatarray().end()) ||
			!_mesh.nodes->IDs->datatarray().back() != _mesh.nodes->IDs->datatarray().size()) {

		reindexElementNodes();
	}
	treindex.end(); timing.addEvent(treindex);
	ESINFO(PROGRESS2) << "Sequential loader:: elements nodes reindexed.";

	timing.totalTime.endWithBarrier();
	timing.printStatsMPI();
}

void SequentialInput::prepareRegions()
{
	size_t threads = environment->OMP_NUM_THREADS;

	for (auto nregion = _meshData.nregions.begin(); nregion != _meshData.nregions.end(); ++nregion) {
		std::vector<size_t> distribution = tarray<eslocal>::distribute(threads, nregion->second.size());
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			auto n = nregion->second.begin();
			if (n != nregion->second.end()) {
				auto nit = std::lower_bound(_meshData.nIDs.begin(), _meshData.nIDs.end(), *n);
				for ( ; n != nregion->second.end(); ++n, ++nit) {
					while (*n != *nit) { ++nit; }
					*n = nit - _meshData.nIDs.begin();
				}
			}
		}
	}

	for (auto eregion = _meshData.eregions.begin(); eregion != _meshData.eregions.end(); ++eregion) {
		std::vector<size_t> distribution = tarray<eslocal>::distribute(threads, eregion->second.size());
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			auto e = eregion->second.begin();
			if (e != eregion->second.end()) {
				auto eit = std::lower_bound(_meshData.eIDs.begin(), _meshData.eIDs.begin(), *e);
				for ( ; e != eregion->second.end(); ++e, ++eit) {
					while (*e != *eit) { ++eit; }
					*e = eit - _meshData.eIDs.begin();
				}
			}
		}
	}
}




