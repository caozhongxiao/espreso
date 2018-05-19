
#include "sequentialinput.h"

#include "../basis/containers/serializededata.h"
#include "../basis/utilities/utils.h"
#include "../basis/utilities/communication.h"
#include "../basis/logging/timeeval.h"

#include "../config/ecf/root.h"

#include "../mesh/mesh.h"
#include "../mesh/preprocessing/meshpreprocessing.h"
#include "../mesh/elements/element.h"
#include "../mesh/store/nodestore.h"
#include "../mesh/store/elementstore.h"
#include "../mesh/store/elementsregionstore.h"
#include "../mesh/store/boundaryregionstore.h"
#include "../old/input/loader.h"

#include <numeric>
#include <algorithm>
#include <fstream>

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

	TimeEvent tnodes("fill nodes"); tnodes.start();
	_nDistribution = { 0, _meshData.nIDs.size() };
	fillNodes();
	tnodes.end(); timing.addEvent(tnodes);
	ESINFO(PROGRESS2) << "Sequential loader:: nodes filled.";

	TimeEvent telements("fill elements"); telements.start();
	_eDistribution = { 0, _meshData.esize.size() };
	fillElements();
	telements.end(); timing.addEvent(telements);
	ESINFO(PROGRESS2) << "Sequential loader:: elements filled.";

	timing.totalTime.endWithBarrier();
	timing.printStatsMPI();
}

void SequentialInput::fillNodes()
{
	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<size_t> cdistribution = tarray<Point>::distribute(threads, _meshData.coordinates.size());

	_mesh.nodes->size = _meshData.coordinates.size();
	_mesh.nodes->distribution = cdistribution;
	_mesh.nodes->IDs = new serializededata<eslocal, eslocal>(1, tarray<eslocal>(cdistribution, _meshData.nIDs));
	_mesh.nodes->coordinates = new serializededata<eslocal, Point>(1, tarray<Point>(cdistribution, _meshData.coordinates));
	_mesh.nodes->ranks = new serializededata<eslocal, int>(1, tarray<int>(threads, _meshData.nIDs.size()));

	_mesh.boundaryRegions.push_back(new BoundaryRegionStore("ALL_NODES", _mesh._eclasses));
	_mesh.boundaryRegions.back()->nodes = new serializededata<eslocal, eslocal>(1, tarray<eslocal>(cdistribution, _meshData.nIDs));

	_mesh.neighboursWithMe.push_back(environment->MPIrank);
}




