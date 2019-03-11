
#include "sequentialinput.h"

#include "basis/containers/serializededata.h"
#include "esinfo/eslog.h"

#include "mesh/mesh.h"
#include "mesh/store/nodestore.h"

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
	eslog::startln("MESIO: BUILD MESH SEQUENTIALLY", "MESIO");

	_nDistribution = { 0, (esint)_meshData.nIDs.size() };
	sortNodes();
	eslog::checkpointln("MESIO: NODES SORTED");

	_eDistribution = { 0, (esint)_meshData.eIDs.size() };
	sortElements();
	eslog::checkpointln("MESIO: ELEMENTS SORTED");

	reindexRegions(_meshData.eregions, _meshData.eIDs);
	reindexRegions(_meshData.nregions, _meshData.nIDs);
	eslog::checkpointln("MESIO: REGION REINDEXED");

	_meshData._nranks.resize(_meshData.nIDs.size());
	_meshData._nrankdist.resize(_meshData.nIDs.size() + 1);
	std::iota(_meshData._nrankdist.begin(), _meshData._nrankdist.end(), 0);
	_mesh.neighboursWithMe.push_back(0);
	eslog::checkpointln("MESIO: RANKS FILLED");

	fillNodes();
	eslog::checkpointln("MESIO: NODES FILLED");

	fillElements();
	eslog::checkpointln("MESIO: ELEMENTS FILLED");

	fillElementRegions();
	fillBoundaryRegions();
	fillNodeRegions();
	eslog::checkpointln("MESIO: REGIONS FILLED");

	if (!(_mesh.nodes->IDs->datatarray().back() != (esint)_mesh.nodes->IDs->datatarray().size())) {
		reindexElementNodes();
		reindexBoundaryNodes();
	}
	eslog::endln("MESIO: ELEMENTS NODES REINDEXED");
}




