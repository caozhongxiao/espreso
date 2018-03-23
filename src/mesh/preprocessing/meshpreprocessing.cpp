
#include "meshpreprocessing.h"

#include "../mesh.h"

#include "../store/store.h"
#include "../store/elementstore.h"
#include "../store/nodestore.h"
#include "../store/elementsregionstore.h"
#include "../store/boundaryregionstore.h"
#include "../elements/element.h"

#include "../../basis/containers/point.h"
#include "../../basis/containers/serializededata.h"
#include "../../basis/matrices/denseMatrix.h"
#include "../../basis/utilities/communication.h"
#include "../../basis/utilities/utils.h"
#include "../../basis/utilities/parser.h"
#include "../../basis/logging/logging.h"
#include "../../basis/logging/timeeval.h"

#include "../../config/ecf/environment.h"

#include <algorithm>
#include <numeric>
#include <cstring>
#include "../store/fetidatastore.h"

using namespace espreso;

size_t MeshPreprocessing::level = 0;

MeshPreprocessing::MeshPreprocessing(Mesh *mesh)
: _mesh(mesh), _morphing(NULL), _timeStatistics(new TimeEval("Mesh preprocessing timing"))
{
	_timeStatistics->totalTime.startWithBarrier();
}

MeshPreprocessing::~MeshPreprocessing()
{
	delete _timeStatistics;
	for (auto it = _timeEvents.begin(); it != _timeEvents.end(); ++it) {
		delete it->second;
	}
}

void MeshPreprocessing::finishPreprocessing()
{
	_timeStatistics->totalTime.endWithBarrier();
	_timeStatistics->printStatsMPI();
}

void MeshPreprocessing::start(const std::string &message)
{
	ESINFO(VERBOSITY(level)) << std::string(2 * level, ' ') << "Mesh preprocessing :: " << message << " started.";
	++level;

	TimeEvent *event;
	if (_timeEvents.find(message) != _timeEvents.end()) {
		event = _timeEvents[message];
	} else {
		_timeEvents[message] = event = new TimeEvent(message);
		_timeStatistics->addPointerToEvent(event);
	}

	event->start();
}

void MeshPreprocessing::skip(const std::string &message)
{
	ESINFO(VERBOSITY(level)) << std::string(2 * level, ' ') << "Mesh preprocessing :: " << message << " skipped.";
}

void MeshPreprocessing::finish(const std::string &message)
{
	TimeEvent *event;
	if (_timeEvents.find(message) != _timeEvents.end()) {
		event = _timeEvents[message];
	} else {
		_timeEvents[message] = event = new TimeEvent(message);
		_timeStatistics->addPointerToEvent(event);
	}

	event->endWithBarrier();

	--level;
	ESINFO(VERBOSITY(level)) << std::string(2 * level, ' ') << "Mesh preprocessing :: " << message << " finished.";
}

void MeshPreprocessing::linkNodesAndElements()
{
	start("link nodes and elements");

	if (_mesh->elements == NULL || _mesh->nodes == NULL) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: fill both elements and nodes.";
	}

	size_t threads = environment->OMP_NUM_THREADS;

	// thread x neighbor x vector(from, to)
	std::vector<std::vector<std::vector<std::pair<eslocal, eslocal> > > > sBuffer(threads);
	std::vector<std::vector<std::pair<eslocal, eslocal> > > rBuffer(_mesh->neighbours.size());
	std::vector<std::pair<eslocal, eslocal> > localLinks;

	localLinks.resize(_mesh->elements->nodes->cend()->begin() - _mesh->elements->nodes->cbegin()->begin());
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto tnodes = _mesh->elements->nodes->cbegin(t);
		size_t offset = _mesh->elements->nodes->cbegin(t)->begin() - _mesh->elements->nodes->cbegin()->begin();

		for (size_t e = _mesh->elements->distribution[t]; e < _mesh->elements->distribution[t + 1]; ++e, ++tnodes) {
			for (auto n = tnodes->begin(); n != tnodes->end(); ++n, ++offset) {
				localLinks[offset].first = *n;
				localLinks[offset].second = _mesh->elements->IDs->datatarray()[e];
			}
		}
	}

	Esutils::sortWithInplaceMerge(localLinks, _mesh->elements->nodes->datatarray().distribution());

	std::vector<size_t> tbegin(threads);
	for (size_t t = 1; t < threads; t++) {
		tbegin[t] = std::lower_bound(localLinks.begin() + tbegin[t - 1], localLinks.end(), _mesh->nodes->distribution[t], [] (std::pair<eslocal, eslocal> &p, size_t n) { return p.first < n; }) - localLinks.begin();
	}

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto ranks = _mesh->nodes->ranks->cbegin(t);
		std::vector<std::vector<std::pair<eslocal, eslocal> > > tBuffer(_mesh->neighbours.size());

		auto begin = localLinks.begin() + tbegin[t];
		auto end = begin;

		auto send = [&] (eslocal id) {
			for (auto it = begin; it != end; ++it) {
				it->first = id;
			}
			size_t i = 0;
			for (auto rank = ranks->begin(); rank != ranks->end(); ++rank) {
				if (*rank != environment->MPIrank) {
					while (_mesh->neighbours[i] < *rank) ++i;
					tBuffer[i].insert(tBuffer[i].end(), begin, end);
				}
			}
		};

		for (size_t n = _mesh->nodes->distribution[t]; n + 1 < _mesh->nodes->distribution[t + 1]; ++n, ++ranks) {
			while (end->first == begin->first) ++end;
			send(_mesh->nodes->IDs->datatarray()[n]);
			begin = end;
		}
		if (t + 1 < threads) {
			while (end != localLinks.end() && end->first == begin->first) ++end;
			send(_mesh->nodes->IDs->datatarray()[_mesh->nodes->distribution[t + 1] - 1]);
		} else {
			end = localLinks.end();
			send(_mesh->nodes->IDs->datatarray().back());
		}

		sBuffer[t].swap(tBuffer);
	}

	for (size_t t = 1; t < threads; t++) {
		for (size_t r = 0; r < sBuffer[0].size(); r++) {
			sBuffer[0][r].insert(sBuffer[0][r].end(), sBuffer[t][r].begin(), sBuffer[t][r].end());
		}
	}

	if (!Communication::exchangeUnknownSize(sBuffer[0], rBuffer, _mesh->neighbours)) {
		ESINFO(ERROR) << "ESPRESO internal error: addLinkFromTo - exchangeUnknownSize.";
	}

	std::vector<size_t> boundaries = { 0, localLinks.size() };
	for (size_t r = 0; r < rBuffer.size(); r++) {
		localLinks.insert(localLinks.end(), rBuffer[r].begin(), rBuffer[r].end());
		boundaries.push_back(localLinks.size());
	}

	Esutils::mergeAppendedData(localLinks, boundaries);

	std::vector<std::vector<eslocal> > linksBoundaries(threads);
	std::vector<std::vector<eslocal> > linksData(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		if (_mesh->nodes->distribution[t] != _mesh->nodes->distribution[t + 1]) {
			auto llink = std::lower_bound(localLinks.begin(), localLinks.end(), _mesh->nodes->IDs->datatarray()[_mesh->nodes->distribution[t]], [] (std::pair<eslocal, eslocal> &p, size_t n) { return p.first < n; });
			eslocal current;

			std::vector<eslocal> tBoundaries, tData;
			if (t == 0) {
				tBoundaries.push_back(0);
			}

			for (size_t n = _mesh->nodes->distribution[t]; n + 1 < _mesh->nodes->distribution[t + 1]; ++n) {
				current = llink->first;
				while (current == llink->first) {
					tData.push_back(llink->second);
					++llink;
				}
				tBoundaries.push_back(llink - localLinks.begin());
			}
			if (t + 1 < threads) {
				current = llink->first;
				while (llink != localLinks.end() && current == llink->first) {
					tData.push_back(llink->second);
					++llink;
				}
			} else {
				while (llink != localLinks.end()) {
					tData.push_back(llink->second);
					++llink;
				}
			}
			tBoundaries.push_back(llink - localLinks.begin());

			linksBoundaries[t].swap(tBoundaries);
			linksData[t].swap(tData);
		}
	}

	_mesh->nodes->elements = new serializededata<eslocal, eslocal>(linksBoundaries, linksData);

	finish("link nodes and elements");
}

void MeshPreprocessing::exchangeHalo()
{
	// halo elements are all elements that have some shared node
	start("exchanging halo");

	if (_mesh->nodes->elements == NULL) {
		this->linkNodesAndElements();
	}

	if (_mesh->elements->regions == NULL) {
		fillRegionMask();
	}

	std::vector<eslocal> eDistribution = _mesh->elements->gatherElementsProcDistribution();

	size_t threads = environment->OMP_NUM_THREADS;
	std::vector<std::vector<eslocal> > sBuffer(_mesh->neighbours.size()), rBuffer(_mesh->neighbours.size());

	std::vector<std::vector<std::vector<eslocal> > > hElements(threads);

	// we have to got through all nodes because intervals are not computed yet
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<std::vector<eslocal> > telements(_mesh->neighbours.size());
		auto elinks = _mesh->nodes->elements->cbegin(t);
		size_t i = 0;

		for (auto ranks = _mesh->nodes->ranks->cbegin(t); ranks != _mesh->nodes->ranks->cend(t); ++ranks, ++elinks) {
			auto begin = elinks->begin();
			auto end = elinks->begin();
			if (ranks->size() > 1) {
				i = 0;
				while (begin != elinks->end() && *begin < eDistribution[environment->MPIrank]) ++begin;
				end = begin;
				while (end != elinks->end() && *end < eDistribution[environment->MPIrank + 1]) ++end;
				for (auto rank = ranks->begin(); rank != ranks->end(); ++rank) {
					if (*rank != environment->MPIrank) {
						while (_mesh->neighbours[i] < *rank) ++i;
						telements[i].insert(telements[i].end(), begin, end);
					}
				}
			}
		}
		hElements[t].swap(telements);
	}

	int rsize = _mesh->elements->regionMaskSize;

	std::vector<std::vector<size_t> > tdist(_mesh->neighbours.size());
	for (size_t n = 0; n < _mesh->neighbours.size(); ++n) {
		tdist[n] = { 0, hElements[0][n].size() };
	}
	for (size_t t = 1; t < threads; t++) {
		for (size_t n = 0; n < _mesh->neighbours.size(); ++n) {
			hElements[0][n].insert(hElements[0][n].end(), hElements[t][n].begin(), hElements[t][n].end());
			tdist[n].push_back(hElements[0][n].size());
		}
	}
	for (size_t n = 0; n < _mesh->neighbours.size(); ++n) {
		Esutils::sortWithInplaceMerge(hElements[0][n], tdist[n]);
	}
	#pragma omp parallel for
	for (size_t n = 0; n < _mesh->neighbours.size(); ++n) {
		Esutils::removeDuplicity(hElements[0][n]);
		tdist[n] = tarray<eslocal>::distribute(threads, hElements[0][n].size());
		sBuffer[n].resize((4 + rsize) * hElements[0][n].size());
	}

	eslocal offset = eDistribution[environment->MPIrank];
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		const auto &IDs = _mesh->elements->IDs->datatarray();
		const auto &body = _mesh->elements->body->datatarray();
		const auto &material = _mesh->elements->material->datatarray();
		const auto &code = _mesh->elements->epointers->datatarray();
		const auto &regions = _mesh->elements->regions->datatarray();
		for (size_t n = 0; n < _mesh->neighbours.size(); ++n) {
			for (size_t e = tdist[n][t]; e < tdist[n][t + 1]; e++) {
				sBuffer[n][(4 + rsize) * e + 0] = IDs[hElements[0][n][e] - offset];
				sBuffer[n][(4 + rsize) * e + 1] = body[hElements[0][n][e] - offset];
				sBuffer[n][(4 + rsize) * e + 2] = material[hElements[0][n][e] - offset];
				sBuffer[n][(4 + rsize) * e + 3] = (eslocal)code[hElements[0][n][e] - offset]->code;
				memcpy(sBuffer[n].data() + (4 + rsize) * e + 4, regions.data() + (hElements[0][n][e] - offset) * rsize, sizeof(int) * rsize);
			}
		}
	}

	if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, _mesh->neighbours)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange halo elements.";
	}

	std::vector<std::vector<eslocal> > hid(threads);
	std::vector<std::vector<int> > hbody(threads), hmaterial(threads), hregions(threads);
	std::vector<std::vector<Element*> > hcode(threads);

	for (size_t n = 0; n < rBuffer.size(); ++n) {
		std::vector<size_t> distribution = tarray<eslocal>::distribute(threads, rBuffer[n].size() / (4 + rsize));
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (size_t e = distribution[t]; e < distribution[t + 1]; ++e) {
				hid[t].push_back(rBuffer[n][(4 + rsize) * e + 0]);
				hbody[t].push_back(rBuffer[n][(4 + rsize) * e + 1]);
				hmaterial[t].push_back(rBuffer[n][(4 + rsize) * e + 2]);
				hcode[t].push_back(_mesh->_eclasses[t] + rBuffer[n][(4 + rsize) * e + 3]);
				hregions[t].insert(hregions[t].end(), rBuffer[n].data() + (4 + rsize) * e + 4, rBuffer[n].data() + (4 + rsize) * e + 4 + rsize);
			}
		}
	}

	_mesh->halo->IDs = new serializededata<eslocal, eslocal>(1, hid);
	_mesh->halo->body = new serializededata<eslocal, int>(1, hbody);
	_mesh->halo->material = new serializededata<eslocal, int>(1, hmaterial);
	_mesh->halo->epointers = new serializededata<eslocal, Element*>(1, hcode);
	_mesh->halo->regions = new serializededata<eslocal, int>(rsize, hregions);

	_mesh->halo->size = _mesh->halo->IDs->datatarray().size();
	_mesh->halo->distribution = _mesh->halo->IDs->datatarray().distribution();

	const auto &hIDs = _mesh->halo->IDs->datatarray();
	std::vector<eslocal> permutation(hIDs.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) { return hIDs[i] < hIDs[j]; });
	_mesh->halo->permute(permutation);

	finish("exchanging halo");
}


void MeshPreprocessing::computeElementsNeighbors()
{
	start("computation of elements neighbors");

	if (_mesh->nodes->elements == NULL) {
		this->linkNodesAndElements();
	}

	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<std::vector<eslocal> > dualDistribution(threads);
	std::vector<std::vector<eslocal> > dualData(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		const auto &IDs = _mesh->elements->IDs->datatarray();
		auto nodes = _mesh->elements->nodes->cbegin(t);
		const auto &epointers = _mesh->elements->epointers->datatarray();

		std::vector<eslocal> ndist = { 0 }, ndata, fdata, tdist, tdata;
		if (t == 0) {
			tdist.push_back(0);
		}
		ndist.reserve(21); // hexa 20
		ndata.reserve(300); // coarse estimation
		fdata.reserve(100);

		for (size_t e = _mesh->elements->distribution[t]; e < _mesh->elements->distribution[t + 1]; ++e, ++nodes) {
			ndist.resize(1);
			ndata.clear();
			for (auto n = nodes->begin(); n != nodes->end(); ++n) {
				auto elements = _mesh->nodes->elements->cbegin() + *n;
				ndist.push_back(ndist.back() + elements->size());
				ndata.insert(ndata.end(), elements->begin(), elements->end());
			}

			for (auto face = epointers[e]->faces->begin(); face != epointers[e]->faces->end(); ++face) {
				fdata.clear();
				for (auto n = face->begin(); n != face->end(); ++n) {
					fdata.insert(fdata.end(), ndata.data() + ndist[*n], ndata.data() + ndist[*n + 1]);
				}
				std::sort(fdata.begin(), fdata.end());

				tdata.push_back(-1);
				auto begin = fdata.begin(), end = begin;
				while (begin != fdata.end()) {
					while (end != fdata.end() && *end == *begin) { ++end; }
					if (*begin != IDs[e] && end - begin == face->size()) {
						tdata.back() = *begin;
						break;
					}
					begin = end;
				}
			}
			tdist.push_back(tdata.size());
		}

		dualDistribution[t].swap(tdist);
		dualData[t].swap(tdata);
	}

	Esutils::threadDistributionToFullDistribution(dualDistribution);

	_mesh->elements->neighbors = new serializededata<eslocal, eslocal>(dualDistribution, dualData);

	finish("computation of elements neighbors");
}

void MeshPreprocessing::computeDecomposedDual(bool separateMaterials, bool separateRegions, bool separateEtype)
{
	start("computation of clusters dual graphs");

	if (_mesh->elements->neighbors == NULL) {
		this->computeElementsNeighbors();
	}

	if (separateRegions && _mesh->elements->regions == NULL) {
		this->fillRegionMask();
	}

	size_t threads = environment->OMP_NUM_THREADS;
	eslocal eBegin = _mesh->elements->gatherElementsProcDistribution()[environment->MPIrank];
	eslocal eEnd   = eBegin + _mesh->elements->size;
	int rsize = _mesh->elements->regionMaskSize;

	std::vector<std::vector<eslocal> > dualDistribution(threads);
	std::vector<std::vector<eslocal> > dualData(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		const auto &body = _mesh->elements->body->datatarray();
		const auto &material = _mesh->elements->material->datatarray();
		const auto &regions = _mesh->elements->regions->datatarray();
		const auto &epointer = _mesh->elements->epointers->datatarray();

		auto neighs = _mesh->elements->neighbors->cbegin(t);

		std::vector<eslocal> edata(20), tdist, tdata;
		if (t == 0) {
			tdist.push_back(0);
		}

		auto areNeighbors = [&] (eslocal e1, eslocal e2) {
			return
					body[e1] == body[e2] &&
					(!separateMaterials || material[e1] == material[e2]) &&
					(!separateRegions || memcmp(regions.data() + e1 * rsize, regions.data() + e2 * rsize, sizeof(int) * rsize) == 0) &&
					(!separateEtype || epointer[e1]->type == epointer[e2]->type);
		};

		for (size_t e = _mesh->elements->distribution[t]; e < _mesh->elements->distribution[t + 1]; ++e, ++neighs) {
			edata.assign(neighs->begin(), neighs->end());
			std::sort(edata.begin(), edata.end());

			for (auto n = edata.begin(); n != edata.end(); ++n) {
				if (eBegin <= *n && *n < eEnd && areNeighbors(e, *n - eBegin)) {
					tdata.push_back(*n - eBegin);
				}
			}

			tdist.push_back(tdata.size());
		}
		dualDistribution[t].swap(tdist);
		dualData[t].swap(tdata);
	}

	Esutils::threadDistributionToFullDistribution(dualDistribution);

	_mesh->elements->decomposedDual = new serializededata<eslocal, eslocal>(dualDistribution, dualData);

	finish("computation of clusters dual graphs");
}

void MeshPreprocessing::computeBoundaryNodes(std::vector<eslocal> &externalBoundary, std::vector<eslocal> &internalBoundary)
{
	start("computation of boundary nodes");

	if (_mesh->elements->neighbors == NULL) {
		this->computeElementsNeighbors();
	}

	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<std::vector<eslocal> > external(threads), internal(threads);

	eslocal eoffset = _mesh->elements->gatherElementsProcDistribution()[environment->MPIrank];

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<eslocal> texternal, tinternal;

		auto neighbors = _mesh->elements->neighbors->cbegin(t);
		auto enodes = _mesh->elements->nodes->cbegin(t);
		for (eslocal d = _mesh->elements->domainDistribution[t]; d < _mesh->elements->domainDistribution[t + 1]; d++) {
			eslocal dbegin = _mesh->elements->elementsDistribution[d];
			eslocal dend = _mesh->elements->elementsDistribution[d + 1];

			for (eslocal e = dbegin; e < dend; ++e, ++neighbors, ++enodes) {
				auto epointer = _mesh->elements->epointers->datatarray()[e];
				auto faces = epointer->faces->begin();

				for (size_t n = 0; n < neighbors->size(); ++n, ++faces) {
					if (neighbors->at(n) == -1) {
						for (auto f = faces->begin(); f != faces->end(); ++f) {
							texternal.push_back(enodes->at(*f));
						}
					} else if (neighbors->at(n) < dbegin + eoffset || dend + eoffset <= neighbors->at(n)) {
						for (auto f = faces->begin(); f != faces->end(); ++f) {
							tinternal.push_back(enodes->at(*f));
						}
					}
				}
			}
		}

		internal[t].swap(tinternal);
		external[t].swap(texternal);
	}

	Esutils::sortWithUniqueMerge(internal);
	Esutils::sortWithUniqueMerge(external);

	externalBoundary.swap(external[0]);

	auto n2i = [ & ] (size_t neighbour) {
		return std::lower_bound(_mesh->neighbours.begin(), _mesh->neighbours.end(), neighbour) - _mesh->neighbours.begin();
	};

	// external nodes need to be synchronized
	std::vector<std::vector<eslocal> > sBuffer(_mesh->neighbours.size()), rBuffer(_mesh->neighbours.size());
	std::vector<eslocal> nExternal;

	for (size_t i = 0; i < externalBoundary.size(); i++) {
		auto nrank = _mesh->nodes->ranks->cbegin() + externalBoundary[i];
		for (auto rank = nrank->begin(); rank != nrank->end(); ++rank) {
			if (*rank != environment->MPIrank) {
				sBuffer[n2i(*rank)].push_back(_mesh->nodes->IDs->datatarray()[externalBoundary[i]]);
			}
		}
	}

	for (size_t n = 0; n < _mesh->neighbours.size(); n++) {
		std::sort(sBuffer[n].begin(), sBuffer[n].end());
	}

	if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, _mesh->neighbours)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange external nodes.";
	}

	for (size_t n = 0; n < _mesh->neighbours.size(); n++) {
		nExternal.insert(nExternal.end(), rBuffer[n].begin(), rBuffer[n].end());
	}
	Esutils::sortAndRemoveDuplicity(nExternal);

	for (size_t n = 0; n < _mesh->neighbours.size(); n++) {
		nExternal.resize(std::set_difference(nExternal.begin(), nExternal.end(), sBuffer[n].begin(), sBuffer[n].end(), nExternal.begin()) - nExternal.begin());
	}

	for (size_t n = 0; n < nExternal.size(); n++) {
		std::vector<std::vector<eslocal> > tnExternal(threads);
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			auto it = std::find(_mesh->nodes->IDs->datatarray().cbegin() + _mesh->nodes->distribution[t], _mesh->nodes->IDs->datatarray().cbegin() + _mesh->nodes->distribution[t + 1], nExternal[n]);
			if (it != _mesh->nodes->IDs->datatarray().cbegin() + _mesh->nodes->distribution[t + 1]) {
				tnExternal[t].push_back(it - _mesh->nodes->IDs->datatarray().cbegin());
			}
		}

		for (size_t t = 0; t < threads; t++) {
			externalBoundary.insert(externalBoundary.end(), tnExternal[t].begin(), tnExternal[t].end());
		}
	}
	std::sort(externalBoundary.begin(), externalBoundary.end());

	internalBoundary.resize(internal[0].size());
	internalBoundary.resize(std::set_difference(internal[0].begin(), internal[0].end(), externalBoundary.begin(), externalBoundary.end(), internalBoundary.begin()) - internalBoundary.begin());

	finish("computation of boundary nodes");
}

void MeshPreprocessing::computeRegionArea(BoundaryRegionStore *store)
{
	double A = 0;
	auto nodes = store->elements->cbegin();
	const auto &epointers = store->epointers->datatarray();
	const auto &coordinates = _mesh->nodes->coordinates->datatarray();
	for (size_t e = 0; e < store->elements->structures(); ++e, ++nodes) {

		DenseMatrix coords(nodes->size(), 3), dND(1, 3);

		const std::vector<DenseMatrix> &dN = *epointers[e]->dN;
		const std::vector<double> &weighFactor = *epointers[e]->weighFactor;

		for (size_t n = 0; n < nodes->size(); ++n) {
			coords(n, 0) = coordinates[nodes->at(n)].x;
			coords(n, 1) = coordinates[nodes->at(n)].y;
			coords(n, 2) = coordinates[nodes->at(n)].z;
		}

		if (store->dimension == 1) {
			for (size_t gp = 0; gp < dN.size(); gp++) {
				dND.multiply(dN[gp], coords);
				A += dND.norm() * weighFactor[gp];
			}
		}
		if (store->dimension == 2) {
			for (size_t gp = 0; gp < dN.size(); gp++) {
				dND.multiply(dN[gp], coords);
				Point v2(dND(0, 0), dND(0, 1), dND(0, 2));
				Point v1(dND(1, 0), dND(1, 1), dND(1, 2));
				Point va = Point::cross(v1, v2);
				A += va.norm() * weighFactor[gp];
			}
		}
	}

	MPI_Allreduce(&A, &store->area, 1, MPI_DOUBLE, MPI_SUM, environment->MPICommunicator);
}
