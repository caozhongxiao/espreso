
#include "sequentialinput.h"

#include "../basis/containers/serializededata.h"
#include "../basis/logging/logging.h"
#include "../basis/logging/timeeval.h"

#include "../mesh/mesh.h"
#include "../mesh/store/nodestore.h"

#include <algorithm>
#include <numeric>

using namespace espreso;

void SequentialInput::buildMesh(PlainMeshData &meshData, Mesh &mesh)
{
	SequentialInput(meshData, mesh);
}

SequentialInput::SequentialInput(PlainMeshData &meshData, Mesh &mesh)
: Input(meshData, mesh)
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

	TimeEvent tpost("reindex regions"); tpost.start();
	reindexRegions(_meshData.eregions, _meshData.eIDs);
	reindexRegions(_meshData.nregions, _meshData.nIDs);
	tpost.end(); timing.addEvent(tpost);
	ESINFO(PROGRESS2) << "Sequential loader:: regions reindexed.";

	TimeEvent tranks("fill ranks"); tranks.start();
	_meshData.nranks.resize(_meshData.nIDs.size());
	_meshData.ndist.resize(_meshData.nIDs.size() + 1);
	std::iota(_meshData.ndist.begin(), _meshData.ndist.end(), 0);
	_mesh.neighboursWithMe.push_back(0);
	tranks.end(); timing.addEvent(tranks);
	ESINFO(PROGRESS2) << "Sequential loader:: ranks filled.";

	TimeEvent tnodes("fill nodes"); tnodes.start();
	fillNodes();
	tnodes.end(); timing.addEvent(tnodes);
	ESINFO(PROGRESS2) << "Sequential loader:: nodes filled.";

	TimeEvent telements("fill elements"); telements.start();
	fillElements();
	telements.end(); timing.addEvent(telements);
	ESINFO(PROGRESS2) << "Sequential loader:: elements filled.";

	TimeEvent tregions("fill regions"); tregions.start();
	fillElementRegions();
	fillBoundaryRegions();
	fillNodeRegions();
	tregions.end(); timing.addEvent(tregions);
	ESINFO(PROGRESS2) << "Sequential loader:: elements filled.";

	TimeEvent treindex("reindex elements nodes"); treindex.start();
	if (!_mesh.nodes->IDs->datatarray().back() != _mesh.nodes->IDs->datatarray().size()) {
		reindexElementNodes();
		reindexBoundaryNodes();
	}
	treindex.end(); timing.addEvent(treindex);
	ESINFO(PROGRESS2) << "Sequential loader:: elements nodes reindexed.";

	timing.totalTime.endWithBarrier();
	timing.printStatsMPI();
}




