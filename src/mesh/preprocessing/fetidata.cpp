
#include "meshpreprocessing.h"

#include "../mesh.h"

#include "../store/store.h"
#include "../store/elementstore.h"
#include "../store/nodestore.h"
#include "../store/elementsregionstore.h"
#include "../store/boundaryregionstore.h"
#include "../store/surfacestore.h"
#include "../elements/element.h"

#include "../../basis/containers/point.h"
#include "../../basis/containers/serializededata.h"
#include "../../basis/matrices/denseMatrix.h"
#include "../../basis/utilities/communication.h"
#include "../../basis/utilities/utils.h"
#include "../../basis/utilities/parser.h"
#include "../../basis/logging/logging.h"

#include "../../config/ecf/environment.h"

#include <algorithm>
#include <numeric>
#include <cstring>
#include "../store/fetidatastore.h"

using namespace espreso;

void MeshPreprocessing::computeSharedFaces()
{
	start("computation of shared faces");

	if (_mesh->nodes->elements == NULL) {
		linkNodesAndElements();
	}

	size_t threads = environment->OMP_NUM_THREADS;
	eslocal eoffset = _mesh->elements->IDs->datatarray().front();

	std::vector<std::vector<eslocal> > inodes(threads);
	std::vector<std::vector<eslocal> > inodesDistribution(threads);
	std::vector<std::vector<std::pair<eslocal, eslocal> > > inodesDomains(threads);

	inodesDistribution.front().push_back(0);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<eslocal> neighDomains;
		std::vector<eslocal> elements;
		std::vector<eslocal> dnodes;
		std::vector<std::pair<eslocal, eslocal> > intervals;
		const auto &epointers = _mesh->elements->epointers->datatarray();
		for (eslocal d = _mesh->elements->domainDistribution[t]; d < _mesh->elements->domainDistribution[t + 1]; d++) {
			neighDomains.clear();
			for (size_t i = 0; i < _mesh->nodes->dintervals[d].size(); i++) {
				auto domains = _mesh->nodes->idomains->cbegin() + _mesh->nodes->dintervals[d][i].pindex;
				neighDomains.insert(neighDomains.end(), domains->begin(), domains->end());
			}
			Esutils::sortAndRemoveDuplicity(neighDomains);

			for (auto nd = neighDomains.begin(); nd != neighDomains.end(); ++nd) {
				if (*nd <= _mesh->elements->firstDomain + d || _mesh->elements->firstDomain + _mesh->elements->ndomains <= *nd) {
					continue;
				}

				intervals.clear();
				for (size_t i = 0; i < _mesh->nodes->dintervals[d].size(); i++) {
					auto domains = _mesh->nodes->idomains->cbegin() + _mesh->nodes->dintervals[d][i].pindex;
					if (std::binary_search(domains->begin(), domains->end(), *nd)) {
						intervals.push_back(std::make_pair(_mesh->nodes->dintervals[d][i].begin, _mesh->nodes->dintervals[d][i].end));
					}
				}
				size_t last = 0;
				for (size_t i = 1; i < intervals.size(); i++) {
					if (intervals[last].second == intervals[i].first) {
						intervals[last].second = intervals[i].second;
					} else {
						++last;
						intervals[last] = intervals[i];
					}
				}
				intervals.resize(last + 1);

				elements.clear();
				for (size_t i = 0; i < intervals.size(); i++) {
					auto nelements = _mesh->nodes->elements->cbegin() + intervals[i].first;
					for (eslocal e = intervals[i].first; e < intervals[i].second; ++e, ++nelements) {
						for (auto ne = nelements->begin(); ne != nelements->end(); ++ne) {
							if (_mesh->elements->elementsDistribution[d] <= *ne - eoffset && *ne - eoffset < _mesh->elements->elementsDistribution[d + 1]) {
								elements.push_back(*ne - eoffset);
							}
						}
					}
				}
				Esutils::sortAndRemoveDuplicity(elements);

				dnodes.clear();
				for (size_t e = 0; e < elements.size(); ++e) {
					auto nodes = _mesh->elements->nodes->cbegin() + elements[e];
					for (auto f = epointers[elements[e]]->faces->cbegin(); f != epointers[elements[e]]->faces->cend(); ++f) {
						size_t size = 0;
						for (auto fn = f->begin(); fn != f->end(); ++fn) {
							auto it = std::lower_bound(intervals.begin(), intervals.end(), nodes->at(*fn), [&] (const std::pair<eslocal, eslocal> &interval, eslocal node) {
								return interval.second < node;
							});
							if (it != intervals.end() && it->first <= nodes->at(*fn) && nodes->at(*fn) < it->second) {
								++size;
							}
						}
						if (size == f->size()) {
							for (auto fn = f->begin(); fn != f->end(); ++fn) {
								dnodes.push_back(nodes->at(*fn));
							}
						}
					}
				}
				Esutils::sortAndRemoveDuplicity(dnodes);
				if (dnodes.size()) {
					inodes[t].insert(inodes[t].end(), dnodes.begin(), dnodes.end());
					inodesDistribution[t].push_back(inodes[t].size());
					inodesDomains[t].push_back(std::make_pair(d, *nd - _mesh->elements->firstDomain));
				}
			}
		}
	}
	Esutils::threadDistributionToFullDistribution(inodesDistribution);
	Esutils::mergeThreadedUniqueData(inodesDistribution);
	Esutils::mergeThreadedUniqueData(inodesDomains);

	if (_mesh->FETIData != NULL) {
		delete _mesh->FETIData;
	}
	_mesh->FETIData = new FETIDataStore();
	_mesh->FETIData->interfaceNodes = new serializededata<eslocal, eslocal>(1, inodes);
	_mesh->FETIData->inodesDistribution = inodesDistribution[0];
	_mesh->FETIData->inodesDomains = inodesDomains[0];

	finish("computation of shared faces");
}

void MeshPreprocessing::computeCornerNodes()
{
	start("computation of corner nodes");

	std::vector<eslocal> domainDistribution = { 0 };
	std::vector<eslocal> domainData;
	auto domains = _mesh->nodes->idomains->cbegin();
	for (size_t i = 0; i < _mesh->nodes->pintervals.size(); ++i, ++domains) {
		if (_mesh->nodes->pintervals[i].end - _mesh->nodes->pintervals[i].begin < 3) {
			for (size_t n = _mesh->nodes->pintervals[i].begin; n < _mesh->nodes->pintervals[i].end; n++) {
				_mesh->FETIData->corners.push_back(n);
				for (auto d = domains->begin(); d != domains->end(); ++d) {
					if (_mesh->elements->firstDomain <= *d && *d < _mesh->elements->firstDomain + _mesh->elements->ndomains) {
						domainData.push_back(*d - _mesh->elements->firstDomain);
					}
				}
				domainDistribution.push_back(domainData.size());
			}
		} else {
			if (domains->size() > 1) {
				_mesh->FETIData->corners.push_back(_mesh->nodes->pintervals[i].begin + (_mesh->nodes->pintervals[i].end - _mesh->nodes->pintervals[i].begin) / 2);
				for (auto d = domains->begin(); d != domains->end(); ++d) {
					if (_mesh->elements->firstDomain <= *d && *d < _mesh->elements->firstDomain + _mesh->elements->ndomains) {
						domainData.push_back(*d - _mesh->elements->firstDomain);
					}
				}
				domainDistribution.push_back(domainData.size());
			}
		}
	}

	_mesh->FETIData->cornerDomains = new serializededata<eslocal, eslocal>(domainDistribution, domainData);

	finish("computation of corner nodes");
}

void MeshPreprocessing::computeDomainsSurface()
{
	start("computation domains surface");

	if (_mesh->elements->dual == NULL) {
		this->computeDual();
	}

	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<eslocal> IDBoundaries = _mesh->elements->gatherElementsDistribution();

	std::vector<std::vector<eslocal> > triangles(threads);
	std::vector<std::vector<eslocal> > intervals(threads);

	intervals.front().push_back(0);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<eslocal> common;
		size_t ncommons, counter;
		eslocal eID = _mesh->elements->distribution[t], eoffset = _mesh->elements->IDs->datatarray().front();
		auto dual = _mesh->elements->dual->cbegin(t);
		auto epointer = _mesh->elements->epointers->cbegin(t);
		auto IDpointer = std::lower_bound(IDBoundaries.begin(), IDBoundaries.end(), eID + eoffset + 1) - 1;
		eslocal begine = *IDpointer, ende = *(IDpointer + 1);

		for (auto e = _mesh->elements->nodes->cbegin(t); e != _mesh->elements->nodes->cend(t); ++e, ++dual, ++epointer, ++eID) {
			if (eID + eoffset >= ende) {
				++IDpointer;
				begine = *IDpointer;
				ende = *(IDpointer + 1);
				intervals[t].push_back(triangles[t].size() / 3);
			}
			if (dual->size() < epointer->front()->faces->structures() || dual->front() < begine || dual->back() >= ende) {

				auto facepointer = epointer->front()->facepointers->datatarray().cbegin();
				for (auto face = epointer->front()->faces->cbegin(); face != epointer->front()->faces->cend(); ++face, ++facepointer) {

					common.clear();
					for (auto n = face->begin(); n != face->end(); ++n) {
						auto nelements = _mesh->nodes->elements->cbegin() + (*e)[*n];
						for (auto ne = nelements->begin(); ne != nelements->end(); ++ne) {
							common.push_back(*ne);
						}
					}
					std::sort(common.begin(), common.end());

					ncommons = counter = 0;
					for (size_t i = 1; i < common.size(); i++) {
						if (common[i - 1] == common[i]) {
							++counter;
						} else {
							if (face->size() == counter + 1) {
								if (begine <= common[i - 1] && common[i - 1] < ende) {
									++ncommons;
								}
							}
							counter = 0;
						}
					}
					if (face->size() == counter + 1) {
						if (begine <= common.back() && common.back() < ende) {
							++ncommons;
						}
					}

					if (ncommons == 1) {
						for (auto n = (*facepointer)->triangles->datatarray().cbegin(); n != (*facepointer)->triangles->datatarray().cend(); ++n) {
							triangles[t].push_back((*e)[face->at(*n)]);
						}
					}
				}
			}
		}
		if (eID + eoffset == ende) {
			intervals[t].push_back(triangles[t].size() / 3);
		}
	}

	if (_mesh->domainsSurface == NULL) {
		_mesh->domainsSurface = new SurfaceStore();
	}

	eslocal tsize = 0;
	for (size_t t = 0; t < threads; t++) {
		for (size_t i = 0; i < intervals[t].size(); i++) {
			intervals[t][i] += tsize;
		}
		tsize += triangles[t].size() / 3;
	}
	Esutils::mergeThreadedUniqueData(intervals);
	Esutils::sortAndRemoveDuplicity(intervals[0]);

	_mesh->domainsSurface->edistribution = intervals[0];

	std::vector<size_t> tdistribution = { 0 };
	for (size_t t = 0; t < threads; t++) {
		tdistribution.push_back(_mesh->domainsSurface->edistribution[_mesh->elements->domainDistribution[t + 1]]);
	}

	serializededata<eslocal, eslocal>::balance(3, triangles, &tdistribution);
	_mesh->domainsSurface->triangles = new serializededata<eslocal, eslocal>(3, triangles);

	std::vector<std::vector<Point> > coordinates(threads);
	std::vector<std::vector<eslocal> > cdistribution(threads);
	cdistribution.front().push_back(0);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (eslocal d = _mesh->elements->domainDistribution[t]; d < _mesh->elements->domainDistribution[t + 1]; d++) {
			std::vector<eslocal> nodes(
					_mesh->domainsSurface->triangles->datatarray().begin() + 3 * _mesh->domainsSurface->edistribution[d],
					_mesh->domainsSurface->triangles->datatarray().begin() + 3 * _mesh->domainsSurface->edistribution[d + 1]);
			Esutils::sortAndRemoveDuplicity(nodes);
			coordinates[t].reserve(coordinates[t].size() + nodes.size());
			for (size_t n = 0; n < nodes.size(); n++) {
				coordinates[t].push_back(_mesh->nodes->coordinates->datatarray()[nodes[n]]);
			}
			for (size_t n = 3 * _mesh->domainsSurface->edistribution[d]; n < 3 * _mesh->domainsSurface->edistribution[d + 1]; n++) {
				_mesh->domainsSurface->triangles->datatarray()[n] = std::lower_bound(nodes.begin(), nodes.end(), _mesh->domainsSurface->triangles->datatarray()[n]) - nodes.begin();
			}
			cdistribution[t].push_back(coordinates[t].size());
		}
	}
	Esutils::threadDistributionToFullDistribution(cdistribution);
	Esutils::mergeThreadedUniqueData(cdistribution);

	_mesh->domainsSurface->coordinates = new serializededata<eslocal, Point>(1, coordinates);
	_mesh->domainsSurface->cdistribution = cdistribution[0];

	finish("computation domains surface");
}

