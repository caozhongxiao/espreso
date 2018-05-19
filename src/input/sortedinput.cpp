
#include "sortedinput.h"

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

void SortedInput::buildMesh(const ECFRoot &configuration, PlainMeshData &meshData, Mesh &mesh)
{
	SortedInput(configuration, meshData, mesh);
}

SortedInput::SortedInput(const ECFRoot &configuration, PlainMeshData &meshData, Mesh &mesh)
: Input(configuration, meshData, mesh)
{
	ESINFO(OVERVIEW) << "Build mesh from sorted elements.";
	TimeEval timing("Load distributed mesh");
	timing.totalTime.startWithBarrier();

	TimeEvent tdistribution("distribute mesh across processes"); tdistribution.start();
	balance();
	checkERegions();
	tdistribution.end(); timing.addEvent(tdistribution);
	ESINFO(PROGRESS2) << "Distributed loader:: data balanced.";

	TimeEvent telements("fill elements"); telements.start();
	fillElements();
	telements.end(); timing.addEvent(telements);
	ESINFO(PROGRESS2) << "Distributed loader:: elements filled.";

	TimeEvent tcoordinates("fill coordinates"); tcoordinates.start();
	fillCoordinates();
	tcoordinates.end(); timing.addEvent(tcoordinates);
	ESINFO(PROGRESS2) << "Distributed loader:: coordinates filled.";

	TimeEvent tnregions("fill node regions"); tnregions.start();
	addNodeRegions();
	tnregions.end(); timing.addEvent(tnregions);
	ESINFO(PROGRESS2) << "Distributed loader:: node regions filled.";

	TimeEvent tbregions("fill boundary regions"); tbregions.start();
	addBoundaryRegions();
	tbregions.end(); timing.addEvent(tbregions);
	ESINFO(PROGRESS2) << "Distributed loader:: boundary regions filled.";

	TimeEvent teregions("fill element regions"); teregions.start();
	addElementRegions();
	teregions.end(); timing.addEvent(teregions);
	ESINFO(PROGRESS2) << "Distributed loader:: elements regions filled.";

	timing.totalTime.endWithBarrier();
	timing.printStatsMPI();
}

void SortedInput::checkERegions()
{
	std::vector<MeshERegion> bregions;

	for (size_t r = 0; r < _meshData.eregions.size(); r++) {
		if (_meshData.eregions[r].min < _eDistribution.back() && _eDistribution.back() < _meshData.eregions[r].max) {
			ESINFO(ERROR) << "ESPRESO Workbench parser error: weird element region.";
		}
		if (_meshData.eregions[r].min >= _eDistribution.back()) {
			bregions.push_back(MeshERegion(std::move(_meshData.eregions[r])));
			_meshData.eregions.erase(_meshData.eregions.begin() + r--);
		}
	}

	size_t bsize = 0;
	std::vector<size_t> rsize = { 0 };
	for (size_t i = 0; i < _meshData.bregions.size(); i++) {
		bsize += _meshData.bregions[i].esize.size();
		rsize.push_back(bsize);
	}

	std::vector<size_t> fdistribution = Communication::getDistribution(bsize, MPITools::operations().sizeToOffsetsSize_t);

	size_t origBSize = _meshData.bregions.size();

	for (size_t r = 0; r < bregions.size(); r++) {
		std::vector<size_t> borders;
		for (int t = 0; t < environment->MPIsize; t++) {
			auto begin = std::lower_bound(bregions[r].elements.begin(), bregions[r].elements.end(), fdistribution[t] + _eDistribution.back());
			auto end = std::lower_bound(bregions[r].elements.begin(), bregions[r].elements.end(), fdistribution[t + 1] + _eDistribution.back());
			if (begin != end) {
				borders.push_back(*begin);
				borders.push_back(borders.back() + end - begin);
			}
		}

		if (!Communication::allGatherUnknownSize(borders)) {
			ESINFO(ERROR) << "ESPRESO internal error: gather bregion borders.";
		}

		bool onlyRename = false;
		for (size_t br = 0; br < origBSize; br++) {
			if (_meshData.bregions[br].min == borders.front() && _meshData.bregions[br].max == borders.back() - 1) {
				_meshData.bregions[br].name = bregions[r].name;
				onlyRename = true;
				break;
			}
		}
		if (onlyRename) {
			continue;
		}

		std::vector<int> tRanks;
		std::vector<std::vector<eslocal> > sBuffer, rBuffer;

		for (int t = 0; t < environment->MPIsize; t++) {
			auto begin = std::lower_bound(bregions[r].elements.begin(), bregions[r].elements.end(), fdistribution[t] + _eDistribution.back());
			auto end = std::lower_bound(bregions[r].elements.begin(), bregions[r].elements.end(), fdistribution[t + 1] + _eDistribution.back());
			if (begin != end) {
				tRanks.push_back(t);
				sBuffer.push_back(std::vector<eslocal>(begin, end));
			}
		}

		if (!Communication::sendVariousTargets(sBuffer, rBuffer, tRanks)) {
			ESINFO(ERROR) << "ESPRESO internal error: send boundary region indices.";
		}

		for (size_t i = 1; i < rBuffer.size(); i++) {
			rBuffer[0].insert(rBuffer[0].end(), rBuffer[i].begin(), rBuffer[i].end());
		}

		auto cmp = [] (EData &edata, eslocal id) {
			return edata.id < id;
		};

		_meshData.bregions.push_back(MeshBRegion());
		_meshData.bregions.back().name = bregions[r].name;
		if (rBuffer.size() && rBuffer.front().size()) {
			for (size_t nr = 0; nr < origBSize; nr++) {
				if (_meshData.bregions[nr].esize.size()) {
					auto begin = std::lower_bound(_meshData.bregions[nr].edata.begin(), _meshData.bregions[nr].edata.end(), rBuffer[0].front(), cmp);
					auto end = std::lower_bound(_meshData.bregions[nr].edata.begin(), _meshData.bregions[nr].edata.end(), rBuffer[0].back() + 1, cmp);
					for (size_t i = begin - _meshData.bregions[nr].edata.begin(), nodes = 0; i < end - _meshData.bregions[nr].edata.begin(); nodes += _meshData.bregions[nr].esize[i++]) {
						_meshData.bregions.back().edata.push_back(_meshData.bregions[nr].edata[i]);
						_meshData.bregions.back().enodes.insert(_meshData.bregions.back().enodes.end(), _meshData.bregions[nr].enodes.begin() + nodes, _meshData.bregions[nr].enodes.begin() + nodes + _meshData.bregions[nr].esize[i]);
						_meshData.bregions.back().esize.push_back(_meshData.bregions[nr].esize[i]);
					}
				}
			}
		}
	}
}

void SortedInput::fillElements()
{
	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<std::vector<eslocal> > tedist(threads);
	std::vector<std::vector<eslocal> > tnodes(threads);
	std::vector<std::vector<eslocal> > eIDs(threads), rData(threads);
	std::vector<std::vector<int> > eMat(threads), eBody(threads);
	std::vector<std::vector<Element*> > epointers(threads);

	std::vector<eslocal> edist = { 0 };
	edist.reserve(_meshData.esize.size() + 1);
	for (size_t e = 0; e < _meshData.esize.size(); e++) {
		edist.push_back(edist.back() + _meshData.esize[e]);
	}

	std::vector<size_t> edistribution = tarray<Point>::distribute(threads, _meshData.esize.size());
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		if(t == 0) {
			tedist[t].insert(tedist[t].end(), edist.begin() + edistribution[t], edist.begin() + edistribution[t + 1] + 1);
		} else {
			tedist[t].insert(tedist[t].end(), edist.begin() + edistribution[t] + 1, edist.begin() + edistribution[t + 1] + 1);
		}
		tnodes[t].insert(tnodes[t].end(), _meshData.enodes.begin() + edist[edistribution[t]], _meshData.enodes.begin() + edist[edistribution[t + 1]]);
		eIDs[t].resize(edistribution[t + 1] - edistribution[t]);
		std::iota(eIDs[t].begin(), eIDs[t].end(), _eDistribution[environment->MPIrank] + edistribution[t]);
		epointers[t].reserve(edistribution[t + 1] - edistribution[t]);
		eBody[t].reserve(edistribution[t + 1] - edistribution[t]);

		for (size_t e = edistribution[t]; e < edistribution[t + 1]; ++e) {
			epointers[t].push_back(&_mesh._eclasses[t][_meshData.edata[e].etype]);
			eBody[t].push_back(_meshData.edata[e].body);
			if (eIDs[t][e - edistribution[t]] != _meshData.edata[e].id) {
				ESINFO(ERROR) << "ESPRESO Workbench parser: not implemented ordering of EBLOCK elements IDs.";
			}
		}


		if (_configuration.input == INPUT_FORMAT::WORKBENCH && _configuration.workbench.keep_material_sets) {
			eMat[t].reserve(edistribution[t + 1] - edistribution[t]);
			for (size_t e = edistribution[t]; e < edistribution[t + 1]; ++e) {
				eMat[t].push_back(_meshData.edata[e].material);
			}
		} else {
			eMat[t].resize(edistribution[t + 1] - edistribution[t]);
		}
	}

	_mesh.elements->size = _meshData.esize.size();
	_mesh.elements->distribution = edistribution;
	_mesh.elements->IDs = new serializededata<eslocal, eslocal>(1, eIDs);
	_mesh.elements->nodes = new serializededata<eslocal, eslocal>(tedist, tnodes);
	_mesh.elements->epointers = new serializededata<eslocal, Element*>(1, epointers);
	_mesh.elements->material = new serializededata<eslocal, int>(1, eMat);
	_mesh.elements->body = new serializededata<eslocal, int>(1, eBody);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		rData[t].resize(edistribution[t + 1] - edistribution[t]);
		std::iota(rData[t].begin(), rData[t].end(), edistribution[t]);
	}
	_mesh.elementsRegions.push_back(new ElementsRegionStore("ALL_ELEMENTS"));
	_mesh.elementsRegions.back()->elements = new serializededata<eslocal, eslocal>(1, rData);
}

void SortedInput::fillCoordinates()
{
	size_t threads = environment->OMP_NUM_THREADS;

	if (environment->MPIsize == 1) {
		std::vector<std::vector<Point> > tcoordinates(threads);
		std::vector<std::vector<eslocal> > nIDs(threads), rData(threads);

		std::vector<size_t> cdistribution = tarray<Point>::distribute(threads, _meshData.coordinates.size());
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			tcoordinates[t].insert(tcoordinates[t].end(), _meshData.coordinates.begin() + cdistribution[t], _meshData.coordinates.begin() + cdistribution[t + 1]);
			nIDs[t].insert(nIDs[t].end(), _meshData.nIDs.begin() + cdistribution[t], _meshData.nIDs.begin() + cdistribution[t + 1]);
		}

		_mesh.nodes->size = _meshData.coordinates.size();
		_mesh.nodes->distribution = cdistribution;
		_mesh.nodes->IDs = new serializededata<eslocal, eslocal>(1, nIDs);
		_mesh.nodes->coordinates = new serializededata<eslocal, Point>(1, tcoordinates);
		_mesh.nodes->ranks = new serializededata<eslocal, int>(1, tarray<int>(threads, _nDistribution.back()));

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			rData[t].resize(cdistribution[t + 1] - cdistribution[t]);
			std::iota(rData[t].begin(), rData[t].end(), cdistribution[t]);
		}
		_mesh.boundaryRegions.push_back(new BoundaryRegionStore("ALL_NODES", _mesh._eclasses));
		_mesh.boundaryRegions.back()->nodes = new serializededata<eslocal, eslocal>(1, rData);

		_mesh.neighboursWithMe.push_back(environment->MPIrank);
		return;
	}

	TimeEval timing("FILL COORDINATES");
	timing.totalTime.startWithBarrier();

	TimeEvent e1("FC SORT ENODES"); e1.start();
	std::vector<std::vector<eslocal> > nodes(threads);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<eslocal> tnodes(_mesh.elements->nodes->datatarray().begin(t), _mesh.elements->nodes->datatarray().end(t));
		Esutils::sortAndRemoveDuplicity(tnodes);
		nodes[t].swap(tnodes);
	}
	Esutils::inplaceMerge(nodes);
	Esutils::removeDuplicity(nodes[0]);

	e1.end(); timing.addEvent(e1);

	TimeEvent e2("FC SBUFFER"); e2.start();

	std::vector<std::vector<eslocal> > sBuffer;
	std::vector<int> sRanks;
	std::vector<int> ssize(environment->MPIsize), rsize(environment->MPIsize);

	for (int t = 0; t < environment->MPIsize; t++) {
		auto begin = std::lower_bound(nodes[0].begin(), nodes[0].end(), _nDistribution[t]);
		auto end = std::lower_bound(nodes[0].begin(), nodes[0].end(), _nDistribution[t + 1]);
		if (end - begin) {
			sBuffer.push_back(std::vector<eslocal>(begin, end));
			sRanks.push_back(t);
		}
		ssize[t] = end - begin;
	}

	e2.end(); timing.addEvent(e2);

	////////

	std::vector<eslocal> rrIDs;

	TimeEvent ee2("FCXX GET RBUFFER SIZES"); ee2.start();

	MPI_Alltoall(ssize.data(), 1, MPI_INT, rsize.data(), 1, MPI_INT, environment->MPICommunicator);

	ee2.end(); timing.addEvent(ee2);

	size_t rrsize = 0;
	for (int t = 0; t < environment->MPIsize; t++) {
		rrsize += rsize[t];
	}
	rrIDs.resize(rrsize);

	TimeEvent ee3("FCXX GET DATA"); ee3.start();

	Communication::allToAllV(nodes[0], rrIDs, ssize, rsize);

	ee3.end(); timing.addEvent(ee3);

	////////

	int avgneighs = 0, nneighs = sRanks.size();
	double allavgsize = 0, avgsize = 0;

	for (size_t i = 0; i < sRanks.size(); i++) {
		avgsize += sBuffer[i].size();
	}
	avgsize /= nneighs;

	MPI_Reduce(&nneighs, &avgneighs, 1, MPI_INT, MPI_SUM, 0, environment->MPICommunicator);
	MPI_Reduce(&avgsize, &allavgsize, 1, MPI_DOUBLE, MPI_SUM, 0, environment->MPICommunicator);

	ESINFO(PROGRESS1) << "FC AVGNEIGHS: " << (double)avgneighs / environment->MPIsize << ", AVGSIZE: " << allavgsize / environment->MPIsize;

	MPI_Reduce(&nneighs, &avgneighs, 1, MPI_INT, MPI_MIN, 0, environment->MPICommunicator);
	MPI_Reduce(&avgsize, &allavgsize, 1, MPI_DOUBLE, MPI_MIN, 0, environment->MPICommunicator);

	ESINFO(PROGRESS1) << "FC MINNEIGHS: " << avgneighs << ", MINSIZE: " << allavgsize;

	MPI_Reduce(&nneighs, &avgneighs, 1, MPI_INT, MPI_MAX, 0, environment->MPICommunicator);
	MPI_Reduce(&avgsize, &allavgsize, 1, MPI_DOUBLE, MPI_MAX, 0, environment->MPICommunicator);

	ESINFO(PROGRESS1) << "FC MAXNEIGHS: " << avgneighs << ", MAXSIZE: " << allavgsize;

	TimeEvent e3("FC EXCHANGE"); e3.start();

	if (!Communication::sendVariousTargets(sBuffer, _rankNodeMap, sRanks, _targetRanks)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange neighbors.";
	}

	e3.end(); timing.addEvent(e3);

	{
		size_t rnodesize = 0;
		for (size_t t = 0; t < _targetRanks.size(); t++) {
			rnodesize += _rankNodeMap[t].size();
		}
		if (rnodesize != rrIDs.size()) {
			ESINFO(ERROR) << "INVALID ALL TO ALL EXCHANGE";
		}
	}

	TimeEvent e4("FC COMPUTE BACKED"); e4.start();

	std::vector<size_t> ndistribution = tarray<Point>::distribute(threads, _meshData.coordinates.size());
	std::vector<std::vector<std::vector<eslocal> > > backedData(threads, std::vector<std::vector<eslocal> >(_targetRanks.size()));

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<eslocal> ranks, ranksOffset;
		std::vector<std::vector<eslocal> > tbackedData(_targetRanks.size());
		std::vector<std::vector<eslocal>::const_iterator> rPointer(_targetRanks.size());

		for (size_t r = 0; r < _targetRanks.size(); r++) {
			rPointer[r] = std::lower_bound(_rankNodeMap[r].begin(), _rankNodeMap[r].end(), _nDistribution[environment->MPIrank] + ndistribution[t]);
		}
		for (size_t n = ndistribution[t]; n < ndistribution[t + 1]; ++n) {
			ranks.clear();
			ranksOffset.clear();
			for (size_t r = 0; r < _targetRanks.size(); r++) {
				if (rPointer[r] != _rankNodeMap[r].end() && *rPointer[r] == _nDistribution[environment->MPIrank] + n) {
					ranksOffset.push_back(r);
					ranks.push_back(_targetRanks[r]);
					++rPointer[r];
				}
			}
			for (size_t r = 0; r < ranks.size(); r++) {
				tbackedData[ranksOffset[r]].push_back(ranksOffset.size());
				tbackedData[ranksOffset[r]].insert(tbackedData[ranksOffset[r]].end(), ranks.begin(), ranks.end());
			}
		}

		backedData[t].swap(tbackedData);
	}

	e4.end(); timing.addEvent(e4);

	#pragma omp parallel for
	for (size_t r = 0; r < _targetRanks.size(); r++) {
		for (size_t t = 1; t < threads; t++) {
			backedData[0][r].insert(backedData[0][r].end(), backedData[t][r].begin(), backedData[t][r].end());
		}
	}

	TimeEvent e5("FC COMPUTE BACKED COORDINATES"); e5.start();

	std::vector<std::vector<Point> > backedCoordinates(_targetRanks.size());
	#pragma omp parallel for
	for (size_t r = 0; r < _targetRanks.size(); r++) {
		backedCoordinates[r].resize(_rankNodeMap[r].size());
	}

	for (size_t r = 0; r < _targetRanks.size(); r++) {
		std::vector<size_t> rdistribution = tarray<eslocal>::distribute(threads, _rankNodeMap[r].size());
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (size_t n = rdistribution[t]; n < rdistribution[t + 1]; ++n) {
				backedCoordinates[r][n] = _meshData.coordinates[_rankNodeMap[r][n] - _nDistribution[environment->MPIrank]];
			}
		}
	}

	e5.end(); timing.addEvent(e5);

	TimeEvent e6("FC RETURN BACKED"); e6.start();

	std::vector<std::vector<eslocal> > nodeRanks(sRanks.size()), allnodes(threads);
	std::vector<std::vector<Point> > coordinates(sRanks.size());

	if (!Communication::sendVariousTargets(backedData[0], nodeRanks, _targetRanks)) {
		ESINFO(ERROR) << "ESPRESO internal error: return node ranks.";
	}
	if (!Communication::sendVariousTargets(backedCoordinates, coordinates, _targetRanks)) {
		ESINFO(ERROR) << "ESPRESO internal error: return coordinates.";
	}

	e6.end(); timing.addEvent(e6);

	TimeEvent e7("FC RANK DATA"); e7.start();

	size_t csize = 0;
	for (size_t i = 0; i < coordinates.size(); i++) {
		csize += coordinates[i].size();
	}

	std::vector<size_t> distribution = tarray<Point>::distribute(threads, csize);
	std::vector<std::vector<eslocal> > rankDistribution(sRanks.size());
	std::vector<std::vector<int> > rankData(sRanks.size());

	#pragma omp parallel for
	for (size_t r = 0; r < sRanks.size(); r++) {
		std::vector<eslocal> trankDistribution;
		std::vector<int> trankData;
		if (r == 0) {
			trankDistribution.push_back(0);
		}

		for (size_t n = 0; n < nodeRanks[r].size(); n += nodeRanks[r][n] + 1) {
			trankData.insert(trankData.end(), nodeRanks[r].begin() + n + 1, nodeRanks[r].begin() + n + 1 + nodeRanks[r][n]);
			trankDistribution.push_back(trankData.size());
		}

		rankDistribution[r].swap(trankDistribution);
		rankData[r].swap(trankData);
	}

	Esutils::threadDistributionToFullDistribution(rankDistribution);

	for (size_t i = threads; i < sRanks.size(); i++) {
		coordinates[threads - 1].insert(coordinates[threads - 1].end(), coordinates[i].begin(), coordinates[i].end());
		rankData[threads - 1].insert(rankData[threads - 1].end(), rankData[i].begin(), rankData[i].end());
		rankDistribution[threads - 1].insert(rankDistribution[threads - 1].end(), rankDistribution[i].begin(), rankDistribution[i].end());
	}
	for (size_t i = threads; i < sRanks.size(); i++) {
		sBuffer[threads - 1].insert(sBuffer[threads - 1].end(), sBuffer[i].begin(), sBuffer[i].end());
	}
	coordinates.resize(threads);
	sBuffer.resize(threads);
	rankData.resize(threads);
	rankDistribution.resize(threads);

	serializededata<eslocal, Point>::balance(1, coordinates, &distribution);
	serializededata<eslocal, eslocal>::balance(1, sBuffer, &distribution);
	serializededata<eslocal, int>::balance(rankDistribution, rankData, &distribution);


	e7.end(); timing.addEvent(e7);

	TimeEvent e8("FC BUILD TARRRAY"); e8.start();

	_mesh.nodes->size = distribution.back();
	_mesh.nodes->distribution = distribution;
	_mesh.nodes->IDs = new serializededata<eslocal, eslocal>(1, sBuffer);
	_mesh.nodes->coordinates = new serializededata<eslocal, Point>(1, coordinates);
	_mesh.nodes->ranks = new serializededata<eslocal, int>(rankDistribution, rankData);

	e8.end(); timing.addEvent(e8);

	TimeEvent e9("FC NEIGHBORS"); e9.start();

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		allnodes[t].resize(distribution[t + 1] - distribution[t]);
		std::iota(allnodes[t].begin(), allnodes[t].end(), distribution[t]);
		Esutils::sortAndRemoveDuplicity(rankData[t]);
	}

	_mesh.boundaryRegions.push_back(new BoundaryRegionStore("ALL_NODES", _mesh._eclasses));
	_mesh.boundaryRegions.back()->nodes = new serializededata<eslocal, eslocal>(1, allnodes);

	for (size_t t = 0; t < threads; t++) {
		_mesh.neighboursWithMe.insert(_mesh.neighboursWithMe.end(), rankData[t].begin(), rankData[t].end());
	}
	Esutils::sortAndRemoveDuplicity(_mesh.neighboursWithMe);

	for (size_t n = 0; n < _mesh.neighboursWithMe.size(); n++) {
		if (_mesh.neighboursWithMe[n] != environment->MPIrank) {
			_mesh.neighbours.push_back(_mesh.neighboursWithMe[n]);
		}
	}

	e9.end(); timing.addEvent(e9);

	TimeEvent e10("FC REINDEX"); e10.start();

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (auto n = _mesh.elements->nodes->begin(t)->begin(); n != _mesh.elements->nodes->end(t)->begin(); ++n) {
			*n = std::lower_bound(_mesh.nodes->IDs->datatarray().begin(), _mesh.nodes->IDs->datatarray().end(), *n) - _mesh.nodes->IDs->datatarray().begin();
		}
	}

	e10.end(); timing.addEvent(e10);

	timing.totalTime.endWithBarrier();
	timing.printStatsMPI();
}

void SortedInput::addNodeRegions()
{
	// assume sorted nodes !!
	size_t threads = environment->OMP_NUM_THREADS;

	if (environment->MPIsize == 1) {
		for (size_t i = 0; i < _meshData.nregions.size(); i++) {
			_mesh.boundaryRegions.push_back(new BoundaryRegionStore(_meshData.nregions[i].name, _mesh._eclasses));

			std::vector<size_t> distribution = tarray<eslocal>::distribute(threads, _meshData.nregions[i].nodes.size());
			std::vector<std::vector<eslocal> > tnodes(threads);
			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				tnodes[t].insert(tnodes[t].end(), _meshData.nregions[i].nodes.begin() + distribution[t], _meshData.nregions[i].nodes.begin() + distribution[t + 1]);
			}
			_mesh.boundaryRegions.back()->nodes = new serializededata<eslocal, eslocal>(1, tnodes);

			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				for (auto n = _mesh.boundaryRegions.back()->nodes->begin(t)->begin(); n != _mesh.boundaryRegions.back()->nodes->end(t)->begin(); ++n) {
					*n = std::lower_bound(_mesh.nodes->IDs->datatarray().begin(), _mesh.nodes->IDs->datatarray().end(), *n) - _mesh.nodes->IDs->datatarray().begin();
				}
			}
		}
		return;
	}

	for (size_t i = 0; i < _meshData.nregions.size(); i++) {
		std::sort(_meshData.nregions[i].nodes.begin(), _meshData.nregions[i].nodes.end());

		std::vector<std::vector<eslocal> > sBuffer, rBuffer;
		std::vector<int> sRanks, tRanks;

		for (int t = 0; t < environment->MPIsize; t++) {
			auto begin = std::lower_bound(_meshData.nregions[i].nodes.begin(), _meshData.nregions[i].nodes.end(), _nDistribution[t]);
			auto end = std::lower_bound(_meshData.nregions[i].nodes.begin(), _meshData.nregions[i].nodes.end(), _nDistribution[t + 1]);
			if (end - begin) {
				sBuffer.push_back(std::vector<eslocal>(begin, end));
				sRanks.push_back(t);
			}
		}

		if (!Communication::sendVariousTargets(sBuffer, rBuffer, sRanks)) {
			ESINFO(ERROR) << "ESPRESO internal error: exchange node region.";
		}

		sBuffer.clear();
		sBuffer.resize(_targetRanks.size());
		for (size_t r = 1; r < rBuffer.size(); r++) {
			rBuffer[0].insert(rBuffer[0].end(), rBuffer[r].begin(), rBuffer[r].end());
		}

		if (rBuffer.size()) {
			#pragma omp parallel for
			for (size_t t = 0; t < _targetRanks.size(); t++) {
				sBuffer[t].resize(rBuffer[0].size());
				sBuffer[t].resize(std::set_intersection(_rankNodeMap[t].begin(), _rankNodeMap[t].end(), rBuffer[0].begin(), rBuffer[0].end(), sBuffer[t].begin()) - sBuffer[t].begin());
			}
		}

		for (size_t t = 0; t < _targetRanks.size(); t++) {
			if (sBuffer[t].size()) {
				tRanks.push_back(t);
			}
		}
		for (size_t t = 0; t < tRanks.size(); t++) {
			sBuffer[t].swap(sBuffer[tRanks[t]]);
			tRanks[t] = _targetRanks[tRanks[t]];
		}
		sBuffer.resize(tRanks.size());

		rBuffer.clear();
		if (!Communication::sendVariousTargets(sBuffer, rBuffer, tRanks)) {
			ESINFO(ERROR) << "ESPRESO internal error: exchange node region to targets.";
		}

		for (size_t t = threads; t < rBuffer.size(); t++) {
			rBuffer[threads - 1].insert(rBuffer[threads - 1].end(), rBuffer[t].begin(), rBuffer[t].end());
		}
		rBuffer.resize(threads);
		serializededata<eslocal, eslocal>::balance(1, rBuffer);

		_mesh.boundaryRegions.push_back(new BoundaryRegionStore(_meshData.nregions[i].name, _mesh._eclasses));
		_mesh.boundaryRegions.back()->nodes = new serializededata<eslocal, eslocal>(1, rBuffer);

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (auto n = _mesh.boundaryRegions.back()->nodes->begin(t)->begin(); n != _mesh.boundaryRegions.back()->nodes->end(t)->begin(); ++n) {
				*n = std::lower_bound(_mesh.nodes->IDs->datatarray().begin(), _mesh.nodes->IDs->datatarray().end(), *n) - _mesh.nodes->IDs->datatarray().begin();
			}
		}
	}
}

void SortedInput::addBoundaryRegions()
{
	size_t threads = environment->OMP_NUM_THREADS;

	if (environment->MPIsize == 1) {
		for (size_t i = 0; i < _meshData.bregions.size(); i++) {
			std::vector<eslocal> edist = { 0 };
			edist.reserve(_meshData.bregions[i].esize.size() + 1);
			for (size_t e = 0; e < _meshData.bregions[i].esize.size(); e++) {
				edist.push_back(edist.back() + _meshData.bregions[i].esize[e]);
			}

			std::vector<std::vector<eslocal> > tedist(threads), tnodes(threads);
			std::vector<std::vector<Element*> > epointers(threads);
			std::vector<size_t> edistribution = tarray<Point>::distribute(threads, _meshData.bregions[i].esize.size());

			tedist.front().push_back(0);
			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				for (size_t n = edistribution[t]; n < edistribution[t + 1]; ++n) {
					tnodes[t].insert(tnodes[t].end(), _meshData.bregions[i].enodes.begin() + edist[n], _meshData.bregions[i].enodes.begin() + edist[n + 1]);
					epointers[t].push_back(&_mesh._eclasses[t][_meshData.bregions[i].edata[n].etype]);
					tedist[t].push_back(tnodes[t].size());
				}
			}

			Esutils::threadDistributionToFullDistribution(tedist);

			_mesh.boundaryRegions.push_back(new BoundaryRegionStore(_meshData.bregions[i].name, _mesh._eclasses));
			_mesh.boundaryRegions.back()->distribution = tarray<eslocal>::distribute(threads, epointers.front().size());
			switch (epointers.front().front()->type) {
			case Element::TYPE::PLANE:
				_mesh.boundaryRegions.back()->dimension = 2;
				break;
			case Element::TYPE::LINE:
				_mesh.boundaryRegions.back()->dimension = 1;
				break;
			default:
				ESINFO(ERROR) << "ESPRESO Workbench parser: invalid boundary region type. Have to be 3D plane or 2D line.";
			}
			_mesh.boundaryRegions.back()->elements = new serializededata<eslocal, eslocal>(tedist, tnodes);
			_mesh.boundaryRegions.back()->epointers = new serializededata<eslocal, Element*>(1, epointers);

			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				for (auto n = _mesh.boundaryRegions.back()->elements->begin(t)->begin(); n != _mesh.boundaryRegions.back()->elements->end(t)->begin(); ++n) {
					*n = std::lower_bound(_mesh.nodes->IDs->datatarray().begin(), _mesh.nodes->IDs->datatarray().end(), *n) - _mesh.nodes->IDs->datatarray().begin();
				}
			}
		}
		return;
	}

	TimeEval timing("BOUNDARY REGIONS");
	timing.totalTime.startWithBarrier();

	TimeEvent e1("BR LINK NODES AND ELEMENTS"); e1.start();

	if (_meshData.bregions.size()) {
		_mesh.preprocessing->linkNodesAndElements();
	}

	e1.end(); timing.addEvent(e1);

	std::vector<eslocal> edistribution = _mesh.elements->gatherElementsProcDistribution();

	for (size_t i = 0; i < _meshData.bregions.size(); i++) {

		TimeEvent e2("BR PREPARE"); e2.start();

		std::vector<eslocal> edist = { 0 };
		edist.reserve(_meshData.bregions[i].esize.size() + 1);
		for (size_t e = 0; e < _meshData.bregions[i].esize.size(); e++) {
			edist.push_back(edist.back() + _meshData.bregions[i].esize[e]);
		}

		std::vector<eslocal> permutation(edist.size() - 1);
		std::iota(permutation.begin(), permutation.end(), 0);
		std::sort(permutation.begin(), permutation.end(), [&] (eslocal e1, eslocal e2) {
			return _meshData.bregions[i].enodes[edist[e1]] < _meshData.bregions[i].enodes[edist[e2]];
		});

		e2.end(); timing.addEvent(e2);

		TimeEvent e3("BR SRANKS"); e3.start();

		std::vector<std::vector<eslocal> > sBuffer, rBuffer;
		std::vector<int> sRanks, tRanks;

		for (int t = 0; t < environment->MPIsize; t++) {
			auto begin = std::lower_bound(permutation.begin(), permutation.end(), _nDistribution[t], [&] (eslocal e, eslocal n) { return _meshData.bregions[i].enodes[edist[e]] < n; });
			auto end = std::lower_bound(permutation.begin(), permutation.end(), _nDistribution[t + 1], [&] (eslocal e, eslocal n) { return _meshData.bregions[i].enodes[edist[e]] < n; });
			if (begin != end) {
				sRanks.push_back(t);
			}
		}
		sBuffer.resize(sRanks.size());

		e3.end(); timing.addEvent(e3);

		TimeEvent e4("BR SBUFFER"); e4.start();

		for (size_t r = 0; r < sRanks.size(); r++) {
			auto begin = std::lower_bound(permutation.begin(), permutation.end(), _nDistribution[sRanks[r]], [&] (eslocal e, eslocal n) { return _meshData.bregions[i].enodes[edist[e]] < n; });
			auto end = std::lower_bound(permutation.begin(), permutation.end(), _nDistribution[sRanks[r] + 1], [&] (eslocal e, eslocal n) { return _meshData.bregions[i].enodes[edist[e]] < n; });
			std::vector<size_t> sdistribution = tarray<eslocal>::distribute(threads, end - begin);
			std::vector<std::vector<eslocal> > tsBuffer(threads);

			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				std::vector<eslocal> ttsBuffer;

				for (auto e = begin + sdistribution[t]; e != begin + sdistribution[t + 1]; ++e) {
					ttsBuffer.push_back(_meshData.bregions[i].edata[*e].etype);
					ttsBuffer.push_back(_meshData.bregions[i].esize[*e]);
					for (eslocal n = 0; n < _meshData.bregions[i].esize[*e]; ++n) {
						ttsBuffer.push_back(_meshData.bregions[i].enodes[edist[*e] + n]);
					}
				}

				tsBuffer[t].swap(ttsBuffer);
			}

			sBuffer[r].push_back(0);
			for (size_t t = 0; t < threads; t++) {
				sBuffer[r].push_back(tsBuffer[t].size() + sBuffer[r].back());
			}
			for (size_t t = 0; t < threads; t++) {
				sBuffer[r].insert(sBuffer[r].end(), tsBuffer[t].begin(), tsBuffer[t].end());
			}

		}

		e4.end(); timing.addEvent(e4);

		int avgneighs = 0, nneighs = sRanks.size();
		double allavgsize = 0, avgsize = 0;
		for (size_t j = 0; j < sRanks.size(); j++) {
			avgsize += sBuffer[j].size();
		}
		avgsize /= nneighs;

		MPI_Reduce(&nneighs, &avgneighs, 1, MPI_INT, MPI_SUM, 0, environment->MPICommunicator);
		MPI_Reduce(&avgsize, &allavgsize, 1, MPI_DOUBLE, MPI_SUM, 0, environment->MPICommunicator);

		ESINFO(PROGRESS1) << "BR AVGNEIGHS: " << (double)avgneighs / environment->MPIsize << ", AVGSIZE: " << allavgsize / environment->MPIsize;

		MPI_Reduce(&nneighs, &avgneighs, 1, MPI_INT, MPI_MIN, 0, environment->MPICommunicator);
		MPI_Reduce(&avgsize, &allavgsize, 1, MPI_DOUBLE, MPI_MIN, 0, environment->MPICommunicator);

		ESINFO(PROGRESS1) << "BR MINNEIGHS: " << avgneighs << ", MINSIZE: " << allavgsize;

		MPI_Reduce(&nneighs, &avgneighs, 1, MPI_INT, MPI_MAX, 0, environment->MPICommunicator);
		MPI_Reduce(&avgsize, &allavgsize, 1, MPI_DOUBLE, MPI_MAX, 0, environment->MPICommunicator);

		ESINFO(PROGRESS1) << "BR MAXNEIGHS: " << avgneighs << ", MAXSIZE: " << allavgsize;

		TimeEvent e5("BR EXCHANGE SBUFFER"); e5.start();

		if (!Communication::sendVariousTargets(sBuffer, rBuffer, sRanks)) {
			ESINFO(ERROR) << "ESPRESO internal error: exchange node region.";
		}

		e5.end(); timing.addEvent(e5);


		nneighs = rBuffer.size();
		avgsize = 0;
		for (size_t j = 0; j < rBuffer.size(); j++) {
			avgsize += rBuffer[j].size();
		}
		avgsize /= nneighs;

		MPI_Reduce(&nneighs, &avgneighs, 1, MPI_INT, MPI_SUM, 0, environment->MPICommunicator);
		MPI_Reduce(&avgsize, &allavgsize, 1, MPI_DOUBLE, MPI_SUM, 0, environment->MPICommunicator);

		ESINFO(PROGRESS1) << "AVGNEIGHS: " << (double)avgneighs / environment->MPIsize << ", AVGSIZE: " << allavgsize / environment->MPIsize;

		MPI_Reduce(&nneighs, &avgneighs, 1, MPI_INT, MPI_MIN, 0, environment->MPICommunicator);
		MPI_Reduce(&avgsize, &allavgsize, 1, MPI_DOUBLE, MPI_MIN, 0, environment->MPICommunicator);

		ESINFO(PROGRESS1) << "MINNEIGHS: " << avgneighs << ", MINSIZE: " << allavgsize;

		MPI_Reduce(&nneighs, &avgneighs, 1, MPI_INT, MPI_MAX, 0, environment->MPICommunicator);
		MPI_Reduce(&avgsize, &allavgsize, 1, MPI_DOUBLE, MPI_MAX, 0, environment->MPICommunicator);

		ESINFO(PROGRESS1) << "MAXNEIGHS: " << avgneighs << ", MAXSIZE: " << allavgsize;

		nneighs = _targetRanks.size();

		MPI_Reduce(&nneighs, &avgneighs, 1, MPI_INT, MPI_MAX, 0, environment->MPICommunicator);
		ESINFO(PROGRESS1) << "AVGTARGETS: " << (double)avgneighs / environment->MPIsize;

		TimeEvent e6("BR PROCESS RBUFFER"); e6.start();

		sBuffer.clear();
		sBuffer.resize(_targetRanks.size());

		for (size_t r = 0; r < rBuffer.size(); r++) {
			std::vector<std::vector<std::vector<eslocal> > > tsBuffer(threads, std::vector<std::vector<eslocal> >(_targetRanks.size()));

			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				std::vector<eslocal> nodes;
				std::vector<std::vector<eslocal> > ttsBuffer(_targetRanks.size());

				for (size_t n = rBuffer[r][t] + threads + 1; n < rBuffer[r][t + 1] + threads + 1; n += 2 + rBuffer[r][n + 1]) {
					nodes.clear();
					nodes.insert(nodes.end(), rBuffer[r].begin() + n + 2, rBuffer[r].begin() + n + 2 + rBuffer[r][n + 1]);
					std::sort(nodes.begin(), nodes.end());
					auto nbegin = std::lower_bound(nodes.begin(), nodes.end(), _nDistribution[environment->MPIrank]);
					auto nend = std::lower_bound(nodes.begin(), nodes.end(), _nDistribution[environment->MPIrank + 1]);

					for (size_t tt = 0; tt < _targetRanks.size(); tt++) {
						auto it = _rankNodeMap[tt].begin();
						bool found = true;
						for (auto current = nbegin; found && current != nend; ++current) {
							it = std::lower_bound(it, _rankNodeMap[tt].end(), *current);
							found = it != _rankNodeMap[t].end() && *it == *current;
						}
						if (found) {
							ttsBuffer[tt].insert(ttsBuffer[tt].end(), rBuffer[r].begin() + n, rBuffer[r].begin() + n + 2 + rBuffer[r][n + 1]);
						}
					}
				}

				tsBuffer[t].swap(ttsBuffer);
			}

			for (size_t tt = 0; tt < _targetRanks.size(); tt++) {
				size_t tsize = 0;
				for (size_t t = 0; t < threads; t++) {
					tsize += tsBuffer[t][tt].size();
				}
				if (tsize) {
					sBuffer[tt].push_back(0);
					for (size_t t = 0; t < threads; t++) {
						sBuffer[tt].push_back(sBuffer[tt].back() + tsBuffer[t][tt].size());
					}
				}
				for (size_t t = 0; t < threads; t++) {
					sBuffer[tt].insert(sBuffer[tt].end(), tsBuffer[t][tt].begin(), tsBuffer[t][tt].end());
				}
			}
		}

		e6.end(); timing.addEvent(e6);

		TimeEvent e7("BR SEND DATA TO POTENTIAL OWNERS"); e7.start();

		for (size_t t = 0; t < _targetRanks.size(); t++) {
			if (sBuffer[t].size()) {
				tRanks.push_back(t);
			}
		}
		for (size_t t = 0; t < tRanks.size(); t++) {
			sBuffer[t].swap(sBuffer[tRanks[t]]);
			tRanks[t] = _targetRanks[tRanks[t]];
		}
		sBuffer.resize(tRanks.size());

		rBuffer.clear();
		if (!Communication::sendVariousTargets(sBuffer, rBuffer, tRanks)) {
			ESINFO(ERROR) << "ESPRESO internal error: exchange node region to targets.";
		}

		e7.end(); timing.addEvent(e7);

		TimeEvent e8("BR BUILD FACES"); e8.start();

		std::vector<std::vector<eslocal> > tedist(threads), tnodes(threads);
		std::vector<std::vector<Element*> > epointers(threads);

		for (size_t r = 0; r < rBuffer.size(); r++) {

			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				std::vector<eslocal> ttedist, ttnodes;
				std::vector<Element*> tepointers;
				if (t == 0 && r == 0) {
					ttedist.push_back(0);
				}
				eslocal foffset = 0;
				if (r && tedist[t].size()) {
					foffset = tedist[t].back();
				}

				std::vector<eslocal> nodes;
				std::vector<eslocal> nlinks;
				int counter;
				bool found = true;
				for (eslocal e = rBuffer[r][t] + threads + 1; e < rBuffer[r][t + 1] + threads + 1; e += 2 + rBuffer[r][e + 1]) {
					found = true;
					for (auto n = rBuffer[r].begin() + e + 2; found && n != rBuffer[r].begin() + e + 2 + rBuffer[r][e + 1]; ++n) {
						auto it = std::lower_bound(_mesh.nodes->IDs->datatarray().begin(), _mesh.nodes->IDs->datatarray().end(), *n);
						if (it != _mesh.nodes->IDs->datatarray().end() && *it == *n) {
							*n = it - _mesh.nodes->IDs->datatarray().begin();
						} else {
							found = false;
						}
					}
					if (found) {
						nlinks.clear();
						for (auto n = rBuffer[r].begin() + e + 2; n != rBuffer[r].begin() + e + 2 + rBuffer[r][e + 1]; ++n) {
							auto links = _mesh.nodes->elements->cbegin() + *n;
							nlinks.insert(nlinks.end(), links->begin(), links->end());
						}
						std::sort(nlinks.begin(), nlinks.end());
						counter = 1;
						for (size_t i = 1; i < nlinks.size(); ++i) {
							if (nlinks[i - 1] == nlinks[i]) {
								++counter;
								if (counter == rBuffer[r][e + 1]) {
									if (_eDistribution[environment->MPIrank] <= nlinks[i] && nlinks[i] < _eDistribution[environment->MPIrank + 1]) {
										ttnodes.insert(ttnodes.end(), rBuffer[r].begin() + e + 2, rBuffer[r].begin() + e + 2 + rBuffer[r][e + 1]);
										ttedist.push_back(ttnodes.size() + foffset);
										tepointers.push_back(&_mesh._eclasses[0][rBuffer[r][e]]);
									}
									break;
								}
							} else {
								counter = 1;
							}
						}
					}
				}

				tedist[t].insert(tedist[t].end(), ttedist.begin(), ttedist.end());
				tnodes[t].insert(tnodes[t].end(), ttnodes.begin(), ttnodes.end());
				epointers[t].insert(epointers[t].end(), tepointers.begin(), tepointers.end());
			}
		}

		e8.end(); timing.addEvent(e8);

		TimeEvent e10("BR CREATE ARRAYS"); e10.start();

		Esutils::threadDistributionToFullDistribution(tedist);

		serializededata<eslocal, eslocal>::balance(tedist, tnodes);
		serializededata<eslocal, Element*>::balance(1, epointers);

		#pragma omp parallel for
		for (size_t t = 1; t < threads; t++) {
			for (size_t e = 0; e < epointers[t].size(); e++) {
				epointers[t][e] = &_mesh._eclasses[t][epointers[t][e] - _mesh._eclasses[0]];
			}
		}

		_mesh.boundaryRegions.push_back(new BoundaryRegionStore(_meshData.bregions[i].name, _mesh._eclasses));
		_mesh.boundaryRegions.back()->dimension = 2;
		if (epointers.front().size()) {
			switch (epointers.front().front()->type) {
			case Element::TYPE::PLANE:
				_mesh.boundaryRegions.back()->dimension = 2;
				break;
			case Element::TYPE::LINE:
				_mesh.boundaryRegions.back()->dimension = 1;
				break;
			default:
				ESINFO(ERROR) << "ESPRESO Workbench parser: invalid boundary region type. Have to be 3D plane or 2D line.";
			}
		}
		int dim = _mesh.boundaryRegions.back()->dimension;
		MPI_Allreduce(&dim, &_mesh.boundaryRegions.back()->dimension, 1, MPI_INT, MPI_MIN, environment->MPICommunicator);

		_mesh.boundaryRegions.back()->elements = new serializededata<eslocal, eslocal>(tedist, tnodes);
		_mesh.boundaryRegions.back()->epointers = new serializededata<eslocal, Element*>(1, epointers);
		_mesh.boundaryRegions.back()->distribution = _mesh.boundaryRegions.back()->epointers->datatarray().distribution();

		e10.end(); timing.addEvent(e10);

		TimeEvent e11("-------"); e11.start();
		e11.end(); timing.addEvent(e11);
	}

	timing.totalTime.endWithBarrier();
	timing.printStatsMPI();
}

void SortedInput::addElementRegions()
{
	size_t threads = environment->OMP_NUM_THREADS;

	if (environment->MPIsize == 1) {
		for (size_t i = 0; i < _meshData.eregions.size(); i++) {
			_mesh.elementsRegions.push_back(new ElementsRegionStore(_meshData.eregions[i].name));

			std::vector<size_t> distribution = tarray<eslocal>::distribute(threads, _meshData.eregions[i].elements.size());
			std::vector<std::vector<eslocal> > telements(threads);
			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				telements[t].insert(telements[t].end(), _meshData.eregions[i].elements.begin() + distribution[t], _meshData.eregions[i].elements.begin() + distribution[t + 1]);
			}
			_mesh.elementsRegions.back()->elements = new serializededata<eslocal, eslocal>(1, telements);
		}
		return;
	}

	for (size_t i = 0; i < _meshData.eregions.size(); i++) {
		std::vector<std::vector<eslocal> > sBuffer, rBuffer;
		std::vector<int> sRanks, tRanks;

		for (int t = 0; t < environment->MPIsize; t++) {
			auto begin = std::lower_bound(_meshData.eregions[i].elements.begin(), _meshData.eregions[i].elements.end(), _eDistribution[t]);
			auto end = std::lower_bound(_meshData.eregions[i].elements.begin(), _meshData.eregions[i].elements.end(), _eDistribution[t + 1]);
			if (end - begin) {
				sBuffer.push_back(std::vector<eslocal>(begin, end));
				sRanks.push_back(t);
			}
		}

		if (!Communication::sendVariousTargets(sBuffer, rBuffer, sRanks)) {
			ESINFO(ERROR) << "ESPRESO internal error: exchange node region.";
		}

		for (size_t t = threads; t < rBuffer.size(); t++) {
			rBuffer[threads - 1].insert(rBuffer[threads - 1].end(), rBuffer[t].begin(), rBuffer[t].end());
		}
		rBuffer.resize(threads);
		serializededata<eslocal, eslocal>::balance(1, rBuffer);

		_mesh.elementsRegions.push_back(new ElementsRegionStore(_meshData.eregions[i].name));
		_mesh.elementsRegions.back()->elements = new serializededata<eslocal, eslocal>(1, rBuffer);

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (auto e = _mesh.elementsRegions.back()->elements->begin(t)->begin(); e != _mesh.elementsRegions.back()->elements->end(t)->begin(); ++e) {
				*e -= _eDistribution[environment->MPIrank];
			}
		}
	}
}
