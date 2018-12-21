
#include "input.h"
#include "workbench/workbench.h"
#include "openfoam/openfoam.h"
#include "meshgenerator/meshgenerator.h"

#include "../basis/containers/serializededata.h"
#include "../basis/logging/logging.h"
#include "../basis/logging/timeeval.h"
#include "../basis/utilities/communication.h"
#include "../basis/utilities/utils.h"

#include "../config/ecf/root.h"
#include "../mesh/mesh.h"
#include "../mesh/elements/element.h"
#include "../mesh/store/nodestore.h"
#include "../mesh/store/elementstore.h"
#include "../mesh/store/elementsregionstore.h"

#include <algorithm>
#include <numeric>

using namespace espreso;

void Input::load(const ECFRoot &configuration, Mesh &mesh)
{
	switch (configuration.input) {
	case INPUT_FORMAT::WORKBENCH:
		WorkbenchLoader::load(configuration.workbench, mesh);
		mesh.update();
		break;
	case INPUT_FORMAT::OPENFOAM:
		OpenFOAMLoader::load(configuration.openfoam, mesh);
		mesh.update();
		break;
	case INPUT_FORMAT::GENERATOR:
	default:
		MeshGenerator::generate(configuration.generator, mesh);
		mesh.update();
		break;
	}
}

void Input::balance()
{
	auto isSorted = [] (const std::vector<esint> &ids) {
		int sorted, allSorted;

		sorted = std::is_sorted(ids.begin(), ids.end());
		MPI_Allreduce(&sorted, &allSorted, 1, MPI_INT, MPI_MIN, environment->MPICommunicator);

		if (!allSorted) {
			return allSorted;
		}

		esint prev = -1, my = ids.size() ? ids.back() : -1;
		if (environment->MPIrank % 2 == 0) {
			if (environment->MPIrank + 1 < environment->MPIsize) {
				MPI_Send(&my, sizeof(esint), MPI_BYTE, environment->MPIrank + 1, 0, environment->MPICommunicator);
			}
			if (environment->MPIrank) {
				MPI_Recv(&prev, sizeof(esint), MPI_BYTE, environment->MPIrank - 1, 0, environment->MPICommunicator, MPI_STATUS_IGNORE);
			}
		} else {
			MPI_Recv(&prev, sizeof(esint), MPI_BYTE, environment->MPIrank - 1, 0, environment->MPICommunicator, MPI_STATUS_IGNORE);
			if (environment->MPIrank + 1 < environment->MPIsize) {
				MPI_Send(&my, sizeof(esint), MPI_BYTE, environment->MPIrank + 1, 0, environment->MPICommunicator);
			}
		}

		if (ids.size()) {
			sorted = prev < ids.front();
		}
		MPI_Allreduce(&sorted, &allSorted, 1, MPI_INT, MPI_MIN, environment->MPICommunicator);
		return allSorted;
	};

	if (isSorted(_meshData.nIDs)) {
		balanceNodes();
	} else {
		balancePermutedNodes();
		sortNodes();
	}

	if (isSorted(_meshData.eIDs)) {
		balanceElements();
	} else {
		balancePermutedElements();
		std::vector<esint> permutation(_meshData.esize.size());
		std::iota(permutation.begin(), permutation.end(), 0);
		std::sort(permutation.begin(), permutation.end(), [&] (esint i, esint j) { return _meshData.eIDs[i] < _meshData.eIDs[j]; });
		sortElements(permutation);
	}
}

std::vector<esint> Input::getDistribution(const std::vector<esint> &IDs, const std::vector<esint> &permutation)
{
	esint myMaxID = 0, maxID;

	if (IDs.size()) {
		myMaxID = IDs[permutation.back()];
	}
	MPI_Allreduce(&myMaxID, &maxID, sizeof(esint), MPI_BYTE, MPITools::esintOperations().max, environment->MPICommunicator);

	std::vector<esint> distribution = tarray<esint>::distribute(environment->MPIsize, maxID + 1);

	return distribution;
}

void Input::balanceNodes()
{
	std::vector<esint> cCurrent = Communication::getDistribution<esint>(_meshData.nIDs.size());
	_nDistribution = tarray<esint>::distribute(environment->MPIsize, cCurrent.back());

	if (!Communication::balance(_meshData.nIDs, cCurrent, _nDistribution)) {
		ESINFO(ERROR) << "ESPRESO internal error: balance node IDs.";
	}
	if (!Communication::balance(_meshData.coordinates, cCurrent, _nDistribution)) {
		ESINFO(ERROR) << "ESPRESO internal error: balance coordinates.";
	}

	esint back = _meshData.nIDs.back();
	MPI_Allgather(&back, sizeof(esint), MPI_BYTE, _nDistribution.data() + 1, sizeof(esint), MPI_BYTE, environment->MPICommunicator);
	for (size_t i = 1; i < _nDistribution.size(); i++) {
		++_nDistribution[i];
	}
}

void Input::balancePermutedNodes()
{
	TimeEval time("BALANCE PERMUTED NODES");
	time.totalTime.startWithBarrier();


	TimeEvent e1("BPN PREPARE DATA");
	e1.start();

	std::vector<esint> permutation(_meshData.nIDs.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (esint i, esint j) { return _meshData.nIDs[i] < _meshData.nIDs[j]; });

	_nDistribution = getDistribution(_meshData.nIDs, permutation);

	std::vector<esint> sBuffer, rBuffer;
	sBuffer.reserve(3 * environment->MPIsize + _meshData.nIDs.size() + _meshData.coordinates.size() * sizeof(Point) / sizeof(esint));

	size_t prevsize;
	auto nbegin = permutation.begin();
	for (int r = 0; r < environment->MPIsize; r++) {
		prevsize = sBuffer.size();
		sBuffer.push_back(0); // total size
		sBuffer.push_back(r); // target
		sBuffer.push_back(0); // number of coordinates

		auto n = nbegin;
		for ( ; n != permutation.end() && _meshData.nIDs[*n] < _nDistribution[r + 1]; ++n) {
			sBuffer.push_back(_meshData.nIDs[*n]);
			sBuffer.insert(sBuffer.end(), reinterpret_cast<const esint*>(_meshData.coordinates.data() + *n), reinterpret_cast<const esint*>(_meshData.coordinates.data() + *n + 1));
		}
		sBuffer[prevsize + 2] = n - nbegin;
		nbegin = n;

		sBuffer[prevsize] = sBuffer.size() - prevsize;
	}

	e1.end();
	time.addEvent(e1);

	TimeEvent e2("BPN EXCHANGE DATA");
	e2.start();

	if (!Communication::allToAllWithDataSizeAndTarget(sBuffer, rBuffer)) {
		ESINFO(ERROR) << "ESPRESO internal error: distribute permuted nodes.";
	}

	e2.end();
	time.addEvent(e2);

	TimeEvent e3("BPN POST PROCESS");
	e3.start();

	_meshData.nIDs.clear();
	_meshData.coordinates.clear();

	size_t offset = 0;
	Point point;
	for (int r = 0; r < environment->MPIsize; r++) {
		++offset;
		size_t csize = rBuffer[++offset]; // coordinates
		++offset;

		for (size_t c = 0; c < csize; ++c) {
			_meshData.nIDs.push_back(rBuffer[offset++]);
			memcpy(&point, rBuffer.data() + offset, sizeof(Point));
			_meshData.coordinates.push_back(point);
			offset += sizeof(Point) / sizeof(esint);
		}
	}

	e3.end();
	time.addEvent(e3);

	time.totalTime.endWithBarrier();
	time.printStatsMPI();
}

void Input::balanceElements()
{
	std::vector<esint> eCurrent = Communication::getDistribution<esint>(_meshData.esize.size());
	_eDistribution = tarray<esint>::distribute(environment->MPIsize, eCurrent.back());

	std::vector<esint> nCurrent = Communication::getDistribution<esint>(_meshData.enodes.size());
	std::vector<esint> nTarget;

	if (environment->MPIrank == 0) {
		nTarget.push_back(0);
	}

	esint nodeOffset = nCurrent[environment->MPIrank];
	size_t eTargetIndex = std::lower_bound(_eDistribution.begin(), _eDistribution.end(), eCurrent[environment->MPIrank] + 1) - _eDistribution.begin();
	for (size_t n = 0; n < _meshData.esize.size(); ++n) {
		nodeOffset += _meshData.esize[n];
		if (eCurrent[environment->MPIrank] + (esint)n + 1 == _eDistribution[eTargetIndex]) {
			nTarget.push_back(nodeOffset);
			++eTargetIndex;
		}
	}
	Communication::allGatherUnknownSize(nTarget);
	nTarget.resize(environment->MPIsize + 1, nTarget.back());

	if (!Communication::balance(_meshData.enodes, nCurrent, nTarget)) {
		ESINFO(ERROR) << "ESPRESO internal error: balance element nodes.";
	}
	if (!Communication::balance(_meshData.esize, eCurrent, _eDistribution)) {
		ESINFO(ERROR) << "ESPRESO internal error: balance element sizes.";
	}
	if (!Communication::balance(_meshData.eIDs, eCurrent, _eDistribution)) {
		ESINFO(ERROR) << "ESPRESO internal error: balance element IDs.";
	}
	if (!Communication::balance(_meshData.body, eCurrent, _eDistribution)) {
		ESINFO(ERROR) << "ESPRESO internal error: balance element bodies.";
	}
	if (!Communication::balance(_meshData.etype, eCurrent, _eDistribution)) {
		ESINFO(ERROR) << "ESPRESO internal error: balance element types.";
	}
	if (!Communication::balance(_meshData.material, eCurrent, _eDistribution)) {
		ESINFO(ERROR) << "ESPRESO internal error: balance element material.";
	}

	auto back = _meshData.eIDs.back();
	MPI_Allgather(&back, sizeof(back), MPI_BYTE, _eDistribution.data() + 1, sizeof(back), MPI_BYTE, environment->MPICommunicator);
	for (size_t i = 1; i < _eDistribution.size(); i++) {
		++_eDistribution[i];
	}
}

void Input::balancePermutedElements()
{
	TimeEval time("BALANCE PERMUTED ELEMENTS");
	time.totalTime.startWithBarrier();

	TimeEvent e1("BPE PREPROCESS");
	e1.start();

	std::vector<esint> permutation(_meshData.eIDs.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (esint i, esint j) { return _meshData.eIDs[i] < _meshData.eIDs[j]; });

	_eDistribution = getDistribution(_meshData.eIDs, permutation);

	std::vector<esint> edist({ 0 });
	edist.reserve(_meshData.esize.size() + 1);
	for (size_t e = 0; e < _meshData.esize.size(); e++) {
		edist.push_back(edist.back() + _meshData.esize[e]);
	}

	std::vector<esint> sBuffer, rBuffer;
	// head, esize, eID, etype, body, material, nodes
	sBuffer.reserve(4 * environment->MPIsize + 5 * _meshData.esize.size() + _meshData.enodes.size());

	size_t prevsize;
	auto ebegin = permutation.begin();
	for (int r = 0; r < environment->MPIsize; r++) {
		prevsize = sBuffer.size();
		sBuffer.push_back(0); // total size
		sBuffer.push_back(r); // target
		sBuffer.push_back(0); // number of elements
		sBuffer.push_back(0); // number of elements nodes

		auto e = ebegin;
		for ( ; e != permutation.end() && _meshData.eIDs[*e] < _eDistribution[r + 1]; ++e) {
			sBuffer.push_back(_meshData.esize[*e]);
			sBuffer.push_back(_meshData.eIDs[*e]);
			sBuffer.push_back(_meshData.etype[*e]);
			sBuffer.push_back(_meshData.body[*e]);
			sBuffer.push_back(_meshData.material[*e]);
			sBuffer.insert(sBuffer.end(), _meshData.enodes.begin() + edist[*e], _meshData.enodes.begin() + edist[*e + 1]);
			sBuffer[prevsize + 3] += edist[*e + 1] - edist[*e];
		}
		sBuffer[prevsize + 2] = e - ebegin;
		ebegin = e;

		sBuffer[prevsize] = sBuffer.size() - prevsize;
	}

	e1.end();
	time.addEvent(e1);

	TimeEvent e2("BPE EXCHANGE");
	e2.start();

	if (!Communication::allToAllWithDataSizeAndTarget(sBuffer, rBuffer)) {
		ESINFO(ERROR) << "ESPRESO internal error: distribute permuted elements.";
	}

	e2.end();
	time.addEvent(e2);

	TimeEvent e3("BPE POST PROCESS");
	e3.start();

	_meshData.esize.clear();
	_meshData.eIDs.clear();
	_meshData.etype.clear();
	_meshData.body.clear();
	_meshData.material.clear();
	_meshData.enodes.clear();

	size_t offset = 0;
	for (int r = 0; r < environment->MPIsize; r++) {
		++offset;
		size_t esize = rBuffer[++offset];
		++offset;
		++offset;

		for (size_t e = 0; e < esize; ++e) {
			_meshData.esize.push_back(rBuffer[offset++]);
			_meshData.eIDs.push_back(rBuffer[offset++]);
			_meshData.etype.push_back(rBuffer[offset++]);
			_meshData.body.push_back(rBuffer[offset++]);
			_meshData.material.push_back(rBuffer[offset++]);
			_meshData.enodes.insert(_meshData.enodes.end(), rBuffer.begin() + offset, rBuffer.begin() + offset + _meshData.esize.back());
			offset += _meshData.esize.back();
		}
	}

	e3.end();
	time.addEvent(e3);

	time.totalTime.endWithBarrier();
	time.printStatsMPI();
}

void Input::sortNodes(bool withElementNodes)
{
	if (std::is_sorted(_meshData.nIDs.begin(), _meshData.nIDs.end())) {
		return;
	}

	std::vector<esint> permutation(_meshData.nIDs.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (esint i, esint j) { return _meshData.nIDs[i] < _meshData.nIDs[j]; });

	std::sort(_meshData.nIDs.begin(), _meshData.nIDs.end());
	Esutils::permute(_meshData.coordinates, permutation);
	Esutils::permute(_nregions, permutation, _nregsize);

	if (_meshData._nranks.size()) {
		std::vector<esint> npermutation(_meshData._nranks.size());
		std::vector<esint> ndist = _meshData._nrankdist;
		for (size_t i = 0, index = 0; i < permutation.size(); i++) {
			_meshData._nrankdist[i + 1] = _meshData._nrankdist[i] + ndist[permutation[i] + 1] - ndist[permutation[i]];
			for (esint n = 0; n < ndist[permutation[i] + 1] - ndist[permutation[i]]; ++n, ++index) {
				npermutation[index] = ndist[permutation[i]] + n;
			}
		}

		Esutils::permute(_meshData._nranks, npermutation);
	}

	if (withElementNodes) {
		std::vector<esint> backpermutation(_meshData.nIDs.size());
		std::iota(backpermutation.begin(), backpermutation.end(), 0);
		std::sort(backpermutation.begin(), backpermutation.end(), [&] (esint i, esint j) { return permutation[i] < permutation[j]; });
		for (size_t n = 0; n < _meshData.enodes.size(); n++) {
			_meshData.enodes[n] = backpermutation[_meshData.enodes[n]];
		}
	}
}

void Input::sortElements()
{
	auto ecomp = [&] (esint i, esint j) {
		if (static_cast<int>(_mesh._eclasses[0][_meshData.etype[i]].type) != static_cast<int>(_mesh._eclasses[0][_meshData.etype[j]].type)) {
			return static_cast<int>(_mesh._eclasses[0][_meshData.etype[i]].type) > static_cast<int>(_mesh._eclasses[0][_meshData.etype[j]].type);
		} else {
			return _meshData.eIDs[i] < _meshData.eIDs[j];
		}
	};

	std::vector<esint> permutation(_meshData.eIDs.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	if (!std::is_sorted(permutation.begin(), permutation.end(), ecomp)) {
		std::sort(permutation.begin(), permutation.end(), ecomp);
		sortElements(permutation);
	}

	for (int type = static_cast<int>(Element::TYPE::VOLUME); type > static_cast<int>(Element::TYPE::POINT); --type) {
		_etypeDistribution.push_back(std::lower_bound(_meshData.etype.begin(), _meshData.etype.end(), type, [&] (int e, int type) {
			return static_cast<int>(_mesh._eclasses[0][e].type) >= type; }) - _meshData.etype.begin()
		);
	}
}

void Input::sortElements(const std::vector<esint> &permutation)
{
	std::vector<esint> edist;
	if (_meshData._edist.size()) {
		edist.swap(_meshData._edist);
	} else {
		edist = std::vector<esint>({ 0 });
		edist.reserve(_meshData.eIDs.size() + 1);
		for (size_t e = 0; e < _meshData.eIDs.size(); e++) {
			edist.push_back(edist.back() + _meshData.esize[e]);
		}
	}

	Esutils::permute(_meshData.eIDs, permutation);
	Esutils::permute(_meshData.esize, permutation);
	Esutils::permute(_meshData.body, permutation);
	Esutils::permute(_meshData.etype, permutation);
	Esutils::permute(_meshData.material, permutation);
	Esutils::permute(_eregions, permutation, _eregsize);

	std::vector<esint> npermutation(_meshData.enodes.size());
	for (size_t i = 0, index = 0; i < permutation.size(); i++) {
		for (esint n = 0; n < _meshData.esize[i]; ++n, ++index) {
			npermutation[index] = edist[permutation[i]] + n;
		}
	}

	Esutils::permute(_meshData.enodes, npermutation);
}

void Input::assignRegions(
		std::map<std::string, std::vector<esint> > &regions, std::vector<esint> &IDs,
		std::vector<esint> &distribution,
		size_t &rsize, std::vector<esint> &rbits)
{
	rsize = regions.size() / (8 * sizeof(esint)) + 1;
	rbits.resize(rsize * IDs.size());

	size_t r = 0;
	for (auto region = regions.begin(); region != regions.end(); ++region, ++r) {
		std::vector<std::vector<esint> > sBuffer, rBuffer;
		std::vector<int> sRanks;

		for (int t = 0; t < environment->MPIsize; t++) {
			auto begin = std::lower_bound(region->second.begin(), region->second.end(), distribution[t]);
			auto end = std::lower_bound(region->second.begin(), region->second.end(), distribution[t + 1]);
			if (end - begin) {
				sBuffer.push_back(std::vector<esint>(begin, end));
				sRanks.push_back(t);
			}
		}

		if (!Communication::sendVariousTargets(sBuffer, rBuffer, sRanks)) {
			ESINFO(ERROR) << "ESPRESO internal error: assign regions.";
		}

		region->second.clear();

		esint byte = r / (8 * sizeof(esint));
		esint bit = 1 << (r % (8 * sizeof(esint)));

		for (size_t t = 0; t < rBuffer.size(); t++) {
			auto it = std::lower_bound(IDs.begin(), IDs.end(), rBuffer[t].front());
			for (size_t i = 0; i < rBuffer[t].size(); ++i) {
				while (rBuffer[t][i] != *it) { ++it; }
				rbits[rsize * (it - IDs.begin()) + byte] |= bit;
			}
		}
	}
}

void Input::fillRegions(std::map<std::string, std::vector<esint> > &regions, size_t &rsize, std::vector<esint> &rbits)
{
	size_t r = 0;
	for (auto region = regions.begin(); region != regions.end(); ++region, ++r) {
		esint byte = r / (8 * sizeof(esint));
		esint bit = 1 << (r % (8 * sizeof(esint)));

		region->second.clear();
		for (size_t i = 0; i < rbits.size() / rsize; ++i) {
			if (rbits[rsize * i + byte] & bit) {
				region->second.push_back(i);
			}
		}
	}
}

void Input::fillNodes()
{
	size_t threads = environment->OMP_NUM_THREADS;

	_mesh.nodes->size = _meshData.coordinates.size();
	_mesh.nodes->distribution = tarray<size_t>::distribute(threads, _meshData.coordinates.size());

	_mesh.nodes->IDs = new serializededata<esint, esint>(1, tarray<esint>(_mesh.nodes->distribution, _meshData.nIDs));
	_mesh.nodes->coordinates = new serializededata<esint, Point>(1, tarray<Point>(_mesh.nodes->distribution, _meshData.coordinates));

	std::vector<size_t> rdistribution = _mesh.nodes->distribution, rdatadistribution = _mesh.nodes->distribution;
	for (size_t t = 1; t < threads; t++) {
		++rdistribution[t];
		rdatadistribution[t] = _meshData._nrankdist[rdistribution[t]];
	}
	++rdistribution[threads];
	rdatadistribution[threads] = _meshData._nrankdist[rdistribution[threads] - 1];

	_mesh.nodes->ranks = new serializededata<esint, int>(tarray<esint>(rdistribution, _meshData._nrankdist), tarray<int>(rdatadistribution, _meshData._nranks));

	_mesh.boundaryRegions.push_back(new BoundaryRegionStore("ALL_NODES", _mesh._eclasses));
	_mesh.boundaryRegions.back()->nodes = new serializededata<esint, esint>(1, tarray<esint>(threads, _meshData.nIDs.size()));
	std::iota(_mesh.boundaryRegions.back()->nodes->datatarray().begin(), _mesh.boundaryRegions.back()->nodes->datatarray().end(), 0);
}

void Input::fillElements()
{
	size_t estart = _mesh.dimension == 3 ? 0 : 1;

	if (!_etypeDistribution.size()) {
		for (int type = static_cast<int>(Element::TYPE::VOLUME); type > static_cast<int>(Element::TYPE::POINT); --type) {
			_etypeDistribution.push_back(std::lower_bound(_meshData.etype.begin(), _meshData.etype.end(), type, [&] (int e, int type) {
				return static_cast<int>(_mesh._eclasses[0][e].type) >= type; }) - _meshData.etype.begin()
			);
		}
	}

	_eDistribution = Communication::getDistribution(_etypeDistribution[estart]);

	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<std::vector<esint> > tedist(threads), tnodes(threads), eIDs(threads), rData(threads);
	std::vector<std::vector<int> > eMat(threads), eBody(threads);
	std::vector<std::vector<Element*> > epointers(threads);

	std::vector<size_t> edistribution = tarray<size_t>::distribute(threads, _etypeDistribution[estart]);

	if (!_meshData._edist.size()) {
		_meshData._edist = { 0 };
		_meshData._edist.reserve(_meshData.esize.size() + 1);
		for (size_t e = 0; e < _meshData.esize.size(); e++) {
			_meshData._edist.push_back(_meshData._edist.back() + _meshData.esize[e]);
		}
	}

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		if(t == 0) {
			tedist[t].insert(tedist[t].end(), _meshData._edist.begin() + edistribution[t], _meshData._edist.begin() + edistribution[t + 1] + 1);
		} else {
			tedist[t].insert(tedist[t].end(), _meshData._edist.begin() + edistribution[t] + 1, _meshData._edist.begin() + edistribution[t + 1] + 1);
		}

		// till now, IDs are irelevant
		eIDs[t].resize(edistribution[t + 1] - edistribution[t]);
		std::iota(eIDs[t].begin(), eIDs[t].end(), _eDistribution[environment->MPIrank] + edistribution[t]);

		eBody[t].insert(eBody[t].end(), _meshData.body.begin() + edistribution[t], _meshData.body.begin() + edistribution[t + 1]);
		tnodes[t].insert(tnodes[t].end(), _meshData.enodes.begin() + _meshData._edist[edistribution[t]], _meshData.enodes.begin() + _meshData._edist[edistribution[t + 1]]);
		eMat[t].insert(eMat[t].end(), _meshData.material.begin() + edistribution[t], _meshData.material.begin() + edistribution[t + 1]);

		epointers[t].resize(edistribution[t + 1] - edistribution[t]);
		for (size_t e = edistribution[t], i = 0; e < edistribution[t + 1]; ++e, ++i) {
			epointers[t][i] = &_mesh._eclasses[t][_meshData.etype[e]];
		}
	}

	_mesh.elements->dimension = _mesh.dimension;
	_mesh.elements->size = _etypeDistribution[estart];
	_mesh.elements->distribution = edistribution;
	_mesh.elements->IDs = new serializededata<esint, esint>(1, eIDs);
	_mesh.elements->procNodes = new serializededata<esint, esint>(tedist, tnodes);
	_mesh.elements->epointers = new serializededata<esint, Element*>(1, epointers);
	_mesh.elements->material = new serializededata<esint, int>(1, eMat);
	_mesh.elements->body = new serializededata<esint, int>(1, eBody);

	_mesh.elementsRegions.push_back(new ElementsRegionStore("ALL_ELEMENTS"));
	_mesh.elementsRegions.back()->elements = new serializededata<esint, esint>(1, tarray<esint>(threads, _mesh.elements->size));
	std::iota(_mesh.elementsRegions.back()->elements->datatarray().begin(), _mesh.elementsRegions.back()->elements->datatarray().end(), 0);
}

void Input::fillNeighbors()
{
	std::vector<int> realnranks = _meshData._nranks;
	Esutils::sortAndRemoveDuplicity(realnranks);

	_mesh.neighboursWithMe.clear();
	_mesh.neighboursWithMe.insert(_mesh.neighboursWithMe.end(), realnranks.begin(), realnranks.end());

	_mesh.neighbours.clear();
	for (size_t n = 0; n < _mesh.neighboursWithMe.size(); n++) {
		if (_mesh.neighboursWithMe[n] != environment->MPIrank) {
			_mesh.neighbours.push_back(_mesh.neighboursWithMe[n]);
		}
	}
}

void Input::fillBoundaryRegions()
{
	size_t threads = environment->OMP_NUM_THREADS;
	size_t estart = _mesh.dimension == 3 ? 0 : 1;

	std::vector<std::vector<esint> > tedist(threads), tnodes(threads);
	std::vector<std::vector<Element*> > epointers(threads);

	if (!_meshData._edist.size()) {
		_meshData._edist = { 0 };
		_meshData._edist.reserve(_meshData.esize.size() + 1);
		for (size_t e = 0; e < _meshData.esize.size(); e++) {
			_meshData._edist.push_back(_meshData._edist.back() + _meshData.esize[e]);
		}
	}

	for (int i = estart; i < 2; i++) {
		std::vector<esint> named;
		for (auto eregion = _meshData.eregions.begin(); eregion != _meshData.eregions.end(); ++eregion) {
			if (eregion->second.size() && _etypeDistribution[estart] <= eregion->second.front()) {
				if (_etypeDistribution[i] <= eregion->second.front() && eregion->second.front() < _etypeDistribution[i + 1]) {
					named.insert(named.end(), eregion->second.begin(), eregion->second.end());
				}
			}
		}
		Esutils::sortAndRemoveDuplicity(named);
		int hasunnamed = 0, add = 0;
		if (named.size() < (size_t)(_etypeDistribution[i + 1] - _etypeDistribution[i])) {
			hasunnamed = 1;
		}
		MPI_Allreduce(&hasunnamed, &add, 1, MPI_INT, MPI_SUM, environment->MPICommunicator);
		if (add) {
			std::vector<esint> &unnamed = i == 1 ? _meshData.eregions["NAMELESS_EDGE_SET"] : _meshData.eregions["NAMELESS_FACE_SET"];
			auto nit = named.begin();
			for (esint index = _etypeDistribution[i]; index < _etypeDistribution[i + 1] && nit != named.end(); ++index) {
				if (*nit != index) {
					unnamed.push_back(index);
				} else {
					++nit;
				}
			}
			size_t prevsize = unnamed.size();
			esint last = named.size() && named.back() > _etypeDistribution[i] ? named.back() + 1 : _etypeDistribution[i];
			unnamed.resize(_etypeDistribution[i + 1] - _etypeDistribution[i] - named.size());
			std::iota(unnamed.begin() + prevsize, unnamed.end(), last);
		}
	}

	for (int i = estart; i < 2; i++) {
		for (auto eregion = _meshData.eregions.begin(); eregion != _meshData.eregions.end(); ++eregion) {
			std::string rname = eregion->first;
			if (StringCompare::caseInsensitivePreffix(_meshData.boundaryprefix, rname)) {
				rname = std::string(rname.begin() + _meshData.boundaryprefix.size(), rname.end());
			}
			int frominterval = 0, add = 0;
			if (eregion->second.size() && _etypeDistribution[i] <= eregion->second.front() && eregion->second.front() < _etypeDistribution[i + 1]) {
				frominterval = 1;
			}
			MPI_Allreduce(&frominterval, &add, 1, MPI_INT, MPI_SUM, environment->MPICommunicator);

			if (add) {
				_mesh.boundaryRegions.push_back(new BoundaryRegionStore(rname, _mesh._eclasses));
				_mesh.boundaryRegions.back()->dimension = 2 - i;

				std::vector<size_t> edistribution = tarray<size_t>::distribute(threads, eregion->second.size());
				std::vector<esint> eregiondist(eregion->second.size() + 1);
				for (size_t e = 0; e < eregion->second.size(); e++) {
					eregiondist[e + 1] = eregiondist[e] + _meshData.esize[eregion->second[e]];
				}

				#pragma omp parallel for
				for (size_t t = 0; t < threads; t++) {
					tedist[t].clear();
					if (t == 0) {
						tedist[t].insert(tedist[t].end(), eregiondist.begin() + edistribution[t], eregiondist.begin() + edistribution[t + 1] + 1);
					} else {
						tedist[t].insert(tedist[t].end(), eregiondist.begin() + edistribution[t] + 1, eregiondist.begin() + edistribution[t + 1] + 1);
					}

					tnodes[t].resize(eregiondist[edistribution[t + 1]] - eregiondist[edistribution[t]]);
					for (size_t e = edistribution[t], index = 0; e < edistribution[t + 1]; ++e) {
						for (esint n = 0; n < _meshData.esize[eregion->second[e]]; ++n, ++index) {
							tnodes[t][index] = _meshData.enodes[_meshData._edist[eregion->second[e]] + n];
						}
					}

					epointers[t].resize(edistribution[t + 1] - edistribution[t]);
					for (size_t e = edistribution[t], i = 0; e < edistribution[t + 1]; ++e, ++i) {
						epointers[t][i] = &_mesh._eclasses[t][_meshData.etype[eregion->second[e]]];
					}
				}

				_mesh.boundaryRegions.back()->distribution = edistribution;
				_mesh.boundaryRegions.back()->procNodes = new serializededata<esint, esint>(tedist, tnodes);
				_mesh.boundaryRegions.back()->epointers = new serializededata<esint, Element*>(1, epointers);
			}
		}
	}
}

void Input::fillNodeRegions()
{
	size_t threads = environment->OMP_NUM_THREADS;

	for (auto nregion = _meshData.nregions.begin(); nregion != _meshData.nregions.end(); ++nregion) {
		_mesh.boundaryRegions.push_back(new BoundaryRegionStore(nregion->first, _mesh._eclasses));
		_mesh.boundaryRegions.back()->nodes = new serializededata<esint, esint>(1, { threads, nregion->second });
	}
}

void Input::fillElementRegions()
{
	size_t threads = environment->OMP_NUM_THREADS;
	size_t estart = _mesh.dimension == 3 ? 0 : 1;

	for (auto eregion = _meshData.eregions.begin(); eregion != _meshData.eregions.end(); ++eregion) {
		int fromelements = 0, add = 0;
		if (eregion->second.size() && eregion->second.front() < _etypeDistribution[estart]) {
			fromelements = 1;
		}
		MPI_Allreduce(&fromelements, &add, 1, MPI_INT, MPI_SUM, environment->MPICommunicator);
		if (add) {
			std::string rname = eregion->first;
			if (StringCompare::caseInsensitivePreffix(_meshData.elementprefix, rname)) {
				rname = std::string(rname.begin() + _meshData.elementprefix.size(), rname.end());
			}
			_mesh.elementsRegions.push_back(new ElementsRegionStore(rname));
			_mesh.elementsRegions.back()->elements = new serializededata<esint, esint>(1, { threads, eregion->second });
		}
	}
}

void Input::reindexRegions(std::map<std::string, std::vector<esint> > &regions, std::vector<esint> &IDs)
{
	for (auto region = regions.begin(); region != regions.end(); ++region) {
		auto n = region->second.begin();
		if (n != region->second.end()) {
			auto nit = std::lower_bound(IDs.begin(), IDs.end(), *n);
			for ( ; n != region->second.end(); ++n, ++nit) {
				while (*n != *nit) { ++nit; }
				*n = nit - IDs.begin();
			}
		}
	}
}

void Input::reindexElementNodes()
{
	size_t threads = environment->OMP_NUM_THREADS;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (auto n = _mesh.elements->procNodes->begin(t)->begin(); n != _mesh.elements->procNodes->end(t)->begin(); ++n) {
			*n = std::lower_bound(_mesh.nodes->IDs->datatarray().begin(), _mesh.nodes->IDs->datatarray().end(), *n) - _mesh.nodes->IDs->datatarray().begin();
		}
	}
}

void Input::reindexBoundaryNodes()
{
	size_t threads = environment->OMP_NUM_THREADS;

	for (size_t r = 0; r < _mesh.boundaryRegions.size(); r++) {
		if (_mesh.boundaryRegions[r]->dimension) {
			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				for (auto n = _mesh.boundaryRegions[r]->procNodes->begin(t)->begin(); n != _mesh.boundaryRegions[r]->procNodes->end(t)->begin(); ++n) {
					*n = std::lower_bound(_mesh.nodes->IDs->datatarray().begin(), _mesh.nodes->IDs->datatarray().end(), *n) - _mesh.nodes->IDs->datatarray().begin();
				}
			}
		}
	}
}



