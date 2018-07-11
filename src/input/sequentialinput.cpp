
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

	TimeEvent tpost("reindex regions"); tpost.start();
	reindexRegions();
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




