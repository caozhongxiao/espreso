
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
#include "../../config/ecf/decomposition.h"

#include "../../wrappers/metis/wmetis.h"
#include "../../wrappers/math/math.h"

#include <algorithm>
#include <numeric>
#include <cstring>
#include "../store/fetidatastore.h"

using namespace espreso;

void MeshPreprocessing::computeSharedFaceNodes()
{
	start("computation of shared face nodes");

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

	if (_mesh->FETIData == NULL) {
		_mesh->FETIData = new FETIDataStore();
	}

	_mesh->FETIData->interfaceNodes = new serializededata<eslocal, eslocal>(1, inodes);
	_mesh->FETIData->inodesDistribution = inodesDistribution[0];
	_mesh->FETIData->inodesDomains = inodesDomains[0];

	finish("computation of shared face nodes");
}

void MeshPreprocessing::computeCornerNodes()
{
	start("computation of corner nodes");

	if (_mesh->FETIData == NULL) {
		_mesh->FETIData = new FETIDataStore();
	}

	std::vector<eslocal> domainDistribution = { 0 };
	std::vector<eslocal> domainData;
	auto domains = _mesh->nodes->idomains->cbegin();
	for (size_t i = 0, ei = 0; i < _mesh->nodes->pintervals.size(); ++i, ++domains) {
		if (domains->size() == 1) {
			continue;
		}
		if (ei < _mesh->nodes->externalIntervals.size() && _mesh->nodes->externalIntervals[ei] == i) {
			++ei;
			eslocal inc = (_mesh->nodes->pintervals[i].end - _mesh->nodes->pintervals[i].begin) / 3;
			inc = inc ? inc : 1;
			for (size_t n = _mesh->nodes->pintervals[i].begin; n < _mesh->nodes->pintervals[i].end; n += inc) {
				_mesh->FETIData->corners.push_back(n);
				for (auto d = domains->begin(); d != domains->end(); ++d) {
					if (_mesh->elements->firstDomain <= *d && *d < _mesh->elements->firstDomain + _mesh->elements->ndomains) {
						domainData.push_back(*d - _mesh->elements->firstDomain);
					}
				}
				domainDistribution.push_back(domainData.size());
			}
		} else {
			_mesh->FETIData->corners.push_back(_mesh->nodes->pintervals[i].begin + (_mesh->nodes->pintervals[i].end - _mesh->nodes->pintervals[i].begin) / 2);
			for (auto d = domains->begin(); d != domains->end(); ++d) {
				if (_mesh->elements->firstDomain <= *d && *d < _mesh->elements->firstDomain + _mesh->elements->ndomains) {
					domainData.push_back(*d - _mesh->elements->firstDomain);
				}
			}
			domainDistribution.push_back(domainData.size());
		}
	}

	_mesh->FETIData->cornerDomains = new serializededata<eslocal, eslocal>(domainDistribution, domainData);

	finish("computation of corner nodes");
}

void MeshPreprocessing::addFixPoints(const serializededata<eslocal, eslocal>* elements, eslocal begin, eslocal end, const serializededata<eslocal, Element*>* epointers, std::vector<eslocal> &fixPoints)
{
	eslocal FIX_POINTS_SIZE = 8;

	auto neighs = [] (std::vector<eslocal> &neighs, Element::CODE code, int node, const eslocal* nodes) {
		switch (code) {
		case Element::CODE::HEXA8:
		case Element::CODE::HEXA20:
			if (node < 4) {
				neighs.push_back((nodes[(node + 1) % 4]));
				neighs.push_back((nodes[(node + 3) % 4]));
				neighs.push_back((nodes[node + 4]));
			} else {
				neighs.push_back((nodes[(node + 1) % 4 + 4]));
				neighs.push_back((nodes[(node + 3) % 4 + 4]));
				neighs.push_back((nodes[node - 4]));
			}
			return 3;
		case Element::CODE::TETRA4:
		case Element::CODE::TETRA10:
			neighs.push_back(nodes[(node + 1) % 4]);
			neighs.push_back(nodes[(node + 2) % 4]);
			neighs.push_back(nodes[(node + 3) % 4]);
			return 3;
		case Element::CODE::PRISMA6:
		case Element::CODE::PRISMA15:
			if (node < 3) {
				neighs.push_back(nodes[(node + 1) % 3]);
				neighs.push_back(nodes[(node + 2) % 3]);
				neighs.push_back(nodes[node + 3]);
			} else {
				neighs.push_back(nodes[(node + 1) % 3 + 3]);
				neighs.push_back(nodes[(node + 2) % 3 + 3]);
				neighs.push_back(nodes[node - 3]);
			}
			return 3;

		case Element::CODE::PYRAMID5:
		case Element::CODE::PYRAMID13:
			if (node == 4) {
				neighs.insert(neighs.end(), nodes, nodes + 4);
				return 4;
			} else {
				neighs.push_back(nodes[(node + 1) % 4]);
				neighs.push_back(nodes[(node + 3) % 4]);
				neighs.push_back(nodes[4]);
				return 3;
			}

		case Element::CODE::TRIANGLE3:
		case Element::CODE::TRIANGLE6:
			neighs.push_back(nodes[(node + 1) % 3]);
			neighs.push_back(nodes[(node + 2) % 3]);
			return 2;

		case Element::CODE::SQUARE4:
		case Element::CODE::SQUARE8:
			neighs.push_back(nodes[(node + 1) % 4]);
			neighs.push_back(nodes[(node + 3) % 4]);
			return 2;

		case Element::CODE::LINE2:
		case Element::CODE::LINE3:
			neighs.push_back(nodes[(node + 1) % 2]);
			return 1;
		}
		return 0;
	};

	std::vector<eslocal> originnodes, neighsnodes;
	auto element = elements->begin() + begin;
	const auto &epointer = epointers->datatarray();
	for (eslocal e = 0; e < end - begin; ++e, ++element) {
		for (int n = 0; n < epointer[begin + e]->coarseNodes; ++n) {
			originnodes.insert(
					originnodes.end(),
					neighs(neighsnodes, epointer[begin + e]->code, n, element->data()),
					element->at(n));
		}
	}

	std::vector<eslocal> permutation(originnodes.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) {
		return originnodes[i] < originnodes[j];
	});

	std::vector<eslocal> ids, dist, data;
	dist.push_back(0);
	ids.push_back(originnodes[permutation[0]]);
	for (size_t i = 0; i < permutation.size(); i++) {
		if (i && originnodes[permutation[i]] != originnodes[permutation[i - 1]]) {
			Esutils::sortAndRemoveDuplicity(data, dist.back());
			dist.push_back(data.size());
			ids.push_back(originnodes[permutation[i]]);
		}
		data.push_back(neighsnodes[permutation[i]]);
	}
	Esutils::sortAndRemoveDuplicity(data, dist.back());
	dist.push_back(data.size());

	for (size_t i = 0; i < data.size(); i++) {
		data[i] = std::lower_bound(ids.begin(), ids.end(), data[i]) - ids.begin();
	}

	std::vector<eslocal> partition(ids.size());
	METISConfiguration options;
	METIS::call(options, ids.size(), dist.data(), data.data(), 0, NULL, NULL, FIX_POINTS_SIZE, partition.data());

	std::vector<std::vector<eslocal> > pids(FIX_POINTS_SIZE), pdist(FIX_POINTS_SIZE, { 0 }), pdata(FIX_POINTS_SIZE);
	for (size_t i = 0; i < partition.size(); i++) {
		pids[partition[i]].push_back(i);
	}
	for (size_t i = 0; i < partition.size(); i++) {
		eslocal p = partition[i];
		for (eslocal j = dist[i]; j < dist[i + 1]; j++) {
			if (partition[data[j]] == p) {
				eslocal index = std::lower_bound(pids[p].begin(), pids[p].end(), data[j]) - pids[p].begin();
				if (pdist[p].size() <= index) {
					pdata[p].push_back(index);
				}
			}
		}
		pdist[p].push_back(pdata[p].size());
	}

	for (eslocal p = 0; p < FIX_POINTS_SIZE; p++) {
		if (pids[p].size()) {
			std::vector<float> vals(pdata[p].size(), 1), x(pids[p].size(), 1. / pids[p].size()), y(pids[p].size());
			float last_l = pids[p].size(), l = 1;

			while (fabs((l - last_l) / l) > 1e-6) {
				MATH::upCSRMatVecProduct(pids[p].size(), pids[p].size(), pdist[p].data(), pdata[p].data(), vals.data(), x.data(), y.data());
				last_l = l;
				l = MATH::vecNorm(pids[p].size(), y.data());
				MATH::vecScale(pids[p].size(), 1 / l, y.data());
				x.swap(y);
			}

			fixPoints.push_back(ids[pids[p][MATH::vecNormMaxIndex(pids[p].size(), x.data())]]);
		}
	}
}

void MeshPreprocessing::computeFixPoints()
{
	start("computation of fix points");

	if (_mesh->FETIData == NULL) {
		_mesh->FETIData = new FETIDataStore();
	}

	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<std::vector<eslocal> > fixPoints(threads), fixPointsDist(threads);

	fixPointsDist.front().push_back(0);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<eslocal> dist, data, partition;
		for (eslocal d = _mesh->elements->domainDistribution[t]; d < _mesh->elements->domainDistribution[t + 1]; d++) {
			size_t size = fixPoints[t].size();
			addFixPoints(_mesh->elements->nodes, _mesh->elements->elementsDistribution[d], _mesh->elements->elementsDistribution[d + 1], _mesh->elements->epointers, fixPoints[t]);
			Esutils::sortAndRemoveDuplicity(fixPoints[t], size);
			fixPointsDist[t].push_back(fixPoints[t].size());
		}
	}

	Esutils::threadDistributionToFullDistribution(fixPointsDist);
	_mesh->FETIData->iFixPointsDistribution.clear();
	for (size_t t = 0; t < threads; t++) {
		_mesh->FETIData->iFixPointsDistribution.insert(_mesh->FETIData->iFixPointsDistribution.end(), fixPointsDist[t].begin(), fixPointsDist[t].end());
		_mesh->FETIData->innerFixPoints.insert(_mesh->FETIData->innerFixPoints.end(), fixPoints[t].begin(), fixPoints[t].end());
	}

	finish("computation of fix points");
}

void MeshPreprocessing::computeFixPointsOnSurface()
{
	start("computation of fix points on surface");

	if (_mesh->domainsSurface == NULL || _mesh->domainsSurface->elements == NULL) {
		computeDomainsSurface();
	}

	if (_mesh->FETIData == NULL) {
		_mesh->FETIData = new FETIDataStore();
	}

	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<std::vector<eslocal> > fixPoints(threads), fixPointsDist(threads);

	fixPointsDist.front().push_back(0);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<eslocal> dist, data, partition;
		for (eslocal d = _mesh->elements->domainDistribution[t]; d < _mesh->elements->domainDistribution[t + 1]; d++) {
			size_t size = fixPoints[t].size();
			addFixPoints(_mesh->domainsSurface->elements, _mesh->domainsSurface->edistribution[d], _mesh->domainsSurface->edistribution[d + 1], _mesh->domainsSurface->epointers, fixPoints[t]);
			Esutils::sortAndRemoveDuplicity(fixPoints[t], size);
			fixPointsDist[t].push_back(fixPoints[t].size());
		}
	}

	Esutils::threadDistributionToFullDistribution(fixPointsDist);
	_mesh->FETIData->sFixPointsDistribution.clear();
	for (size_t t = 0; t < threads; t++) {
		_mesh->FETIData->sFixPointsDistribution.insert(_mesh->FETIData->sFixPointsDistribution.end(), fixPointsDist[t].begin(), fixPointsDist[t].end());
		_mesh->FETIData->surfaceFixPoints.insert(_mesh->FETIData->surfaceFixPoints.end(), fixPoints[t].begin(), fixPoints[t].end());
	}

	finish("computation of fix points on surface");
}

void MeshPreprocessing::computeDomainsSurface()
{
	start("computation domains surface");

	if (_mesh->elements->neighbors == NULL) {
		this->computeElementsNeighbors();
	}

	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<std::vector<eslocal> > faces(threads), facesDistribution(threads), ecounter(threads, std::vector<eslocal>((int)Element::CODE::SIZE));
	std::vector<std::vector<Element*> > fpointer(threads);
	std::vector<std::vector<size_t> > intervals(threads);

	eslocal eoffset = _mesh->elements->gatherElementsProcDistribution()[environment->MPIrank];

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<eslocal> tfaces, tfacesDistribution, tecounter((int)Element::CODE::SIZE);
		std::vector<Element*> tfpointer;
		std::vector<size_t> tintervals;
		if (t == 0) {
			tfacesDistribution.push_back(0);
			tintervals.push_back(0);
		}

		auto neighbors = _mesh->elements->neighbors->cbegin(t);
		auto enodes = _mesh->elements->nodes->cbegin(t);
		for (eslocal d = _mesh->elements->domainDistribution[t]; d < _mesh->elements->domainDistribution[t + 1]; d++) {
			eslocal dbegin = _mesh->elements->elementsDistribution[d];
			eslocal dend = _mesh->elements->elementsDistribution[d + 1];

			for (eslocal e = dbegin; e < dend; ++e, ++neighbors, ++enodes) {
				auto epointer = _mesh->elements->epointers->datatarray()[e];
				auto faces = epointer->faces->begin();
				auto facepointer = epointer->facepointers->datatarray().begin();

				for (size_t n = 0; n < neighbors->size(); ++n, ++faces, ++facepointer) {
					if (neighbors->at(n) < dbegin + eoffset || dend + eoffset <= neighbors->at(n)) {
						for (auto f = faces->begin(); f != faces->end(); ++f) {
							tfaces.push_back(enodes->at(*f));
						}
						tfacesDistribution.push_back(tfaces.size());
						tfpointer.push_back(*facepointer);
						++tecounter[(int)(*facepointer)->code];
					}
				}
			}

			tintervals.push_back(tfacesDistribution.size());
		}

		faces[t].swap(tfaces);
		facesDistribution[t].swap(tfacesDistribution);
		fpointer[t].swap(tfpointer);
		ecounter[t].swap(tecounter);
		intervals[t].swap(tintervals);
	}

	if (_mesh->domainsSurface == NULL) {
		_mesh->domainsSurface = new SurfaceStore();
	}

	for (size_t t = 1; t < threads; t++) {
		for (size_t e = 0; e < ecounter[0].size(); e++) {
			ecounter[0][e] += ecounter[t][e];
		}
	}

	_mesh->domainsSurface->epointers = new serializededata<eslocal, Element*>(1, fpointer);
	_mesh->domainsSurface->ecounters = ecounter[0];

	for (size_t i = 1; i < intervals[0].size(); i++) {
		--intervals[0][i];
	}
	eslocal tsize = facesDistribution[0].size() - 1;
	for (size_t t = 1; t < threads; t++) {
		for (size_t i = 0; i < intervals[t].size(); i++) {
			intervals[t][i] += tsize;
		}
		tsize += facesDistribution[t].size();
	}
	Esutils::mergeThreadedUniqueData(intervals);
	Esutils::sortAndRemoveDuplicity(intervals[0]);

	_mesh->domainsSurface->edistribution = intervals[0];

	std::vector<size_t> tdistribution = { 0 };
	for (size_t t = 0; t < threads; t++) {
		tdistribution.push_back(_mesh->domainsSurface->edistribution[_mesh->elements->domainDistribution[t + 1]]);
	}

	if (ecounter[0][(int)Element::CODE::TRIANGLE3] == _mesh->domainsSurface->edistribution.back()) {
		serializededata<eslocal, eslocal>::balance(3, faces, &tdistribution);
		_mesh->domainsSurface->elements = new serializededata<eslocal, eslocal>(3, faces);
		_mesh->domainsSurface->triangles = _mesh->domainsSurface->elements;
		_mesh->domainsSurface->tdistribution = _mesh->domainsSurface->edistribution;
	} else {
		Esutils::threadDistributionToFullDistribution(facesDistribution);
		serializededata<eslocal, eslocal>::balance(facesDistribution, faces, &tdistribution);
		_mesh->domainsSurface->elements = new serializededata<eslocal, eslocal>(facesDistribution, faces);
	}

	std::vector<std::vector<Point> > coordinates(threads);
	std::vector<std::vector<size_t> > cdistribution(threads);

	std::vector<eslocal> sintervals = _mesh->nodes->externalIntervals;
	auto idomains = _mesh->nodes->idomains->cbegin();
	for (size_t i = 0; i < _mesh->nodes->pintervals.size(); ++i, ++idomains) {
		if (idomains->size() > 1) {
			sintervals.push_back(i);
		}
	}
	Esutils::sortAndRemoveDuplicity(sintervals);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<Point> tcoordinates;
		std::vector<size_t> tcdistribution;
		std::vector<eslocal> ibounds, ioffset;
		if (t == 0) {
			tcdistribution.push_back(0);
		}

		size_t nsize = 0;
		for (eslocal d = _mesh->elements->domainDistribution[t]; d < _mesh->elements->domainDistribution[t + 1]; d++) {
			auto sit = sintervals.begin();
			for (size_t i = 0; i < _mesh->nodes->dintervals[d].size(); ++i) {
				while (*sit < _mesh->nodes->dintervals[d][i].pindex) { ++sit; }
				if (*sit == _mesh->nodes->dintervals[d][i].pindex) {
					nsize += _mesh->nodes->dintervals[d][i].end - _mesh->nodes->dintervals[d][i].begin;
				}
			}
		}
		tcoordinates.reserve(nsize);
		for (eslocal d = _mesh->elements->domainDistribution[t]; d < _mesh->elements->domainDistribution[t + 1]; d++) {
			auto sit = sintervals.begin();
			for (size_t i = 0; i < _mesh->nodes->dintervals[d].size(); ++i) {
				while (*sit < _mesh->nodes->dintervals[d][i].pindex) { ++sit; }
				if (*sit == _mesh->nodes->dintervals[d][i].pindex) {
					tcoordinates.insert(tcoordinates.end(),
							_mesh->nodes->coordinates->datatarray().data() + _mesh->nodes->dintervals[d][i].begin,
							_mesh->nodes->coordinates->datatarray().data() + _mesh->nodes->dintervals[d][i].end);
				}
			}
			tcdistribution.push_back(tcoordinates.size());
		}

		for (eslocal d = _mesh->elements->domainDistribution[t]; d < _mesh->elements->domainDistribution[t + 1]; d++) {
			auto sit = sintervals.begin();
			ibounds.clear();
			ioffset = { 0 };
			for (size_t i = 0; i < _mesh->nodes->dintervals[d].size(); ++i) {
				while (*sit < _mesh->nodes->dintervals[d][i].pindex) { ++sit; }
				if (*sit == _mesh->nodes->dintervals[d][i].pindex) {
					ibounds.push_back(_mesh->nodes->dintervals[d][i].begin);
					ioffset.push_back(ioffset.back() + _mesh->nodes->dintervals[d][i].end - _mesh->nodes->dintervals[d][i].begin);
				}
			}

			for (
					auto n = (_mesh->domainsSurface->elements->begin() + _mesh->domainsSurface->edistribution[d])->begin();
					n != (_mesh->domainsSurface->elements->begin() + _mesh->domainsSurface->edistribution[d + 1])->begin();
					++n) {

				size_t i = std::lower_bound(ibounds.begin(), ibounds.end(), *n + 1) - ibounds.begin() - 1;
				*n = *n - ibounds[i] + ioffset[i];
			}
		}

		cdistribution[t].swap(tcdistribution);
		coordinates[t].swap(tcoordinates);
	}
	Esutils::threadDistributionToFullDistribution(cdistribution);
	Esutils::mergeThreadedUniqueData(cdistribution);

	_mesh->domainsSurface->coordinates = new serializededata<eslocal, Point>(1, coordinates);
	_mesh->domainsSurface->cdistribution = cdistribution[0];

	finish("computation domains surface");
}

void MeshPreprocessing::triangularizeDomainSurface()
{
	start("triangularize domain surface");

	if (_mesh->domainsSurface->elements == NULL) {
		computeDomainsSurface();
	}

	size_t threads = environment->OMP_NUM_THREADS;

	if (_mesh->domainsSurface->triangles == NULL) {

		std::vector<std::vector<eslocal> > triangles(threads);
		std::vector<std::vector<size_t> > intervals(threads);

		intervals.front().push_back(0);
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (eslocal d = _mesh->elements->domainDistribution[t]; d < _mesh->elements->domainDistribution[t + 1]; d++) {
				auto elements = _mesh->domainsSurface->elements->cbegin() + _mesh->domainsSurface->edistribution[d];
				auto epointer = _mesh->domainsSurface->epointers->datatarray().begin() + _mesh->domainsSurface->edistribution[d];

				for (eslocal e = _mesh->domainsSurface->edistribution[d]; e < _mesh->domainsSurface->edistribution[d + 1]; ++e, ++elements, ++epointer) {
					for (auto n = (*epointer)->triangles->datatarray().cbegin(); n != (*epointer)->triangles->datatarray().cend(); ++n) {
						triangles[t].push_back(elements->at(*n));
					}
				}
				intervals[t].push_back(triangles[t].size() / 3);
			}
		}

		Esutils::threadDistributionToFullDistribution(intervals);
		Esutils::mergeThreadedUniqueData(intervals);

		_mesh->domainsSurface->tdistribution = intervals[0];
		_mesh->domainsSurface->triangles = new serializededata<eslocal, eslocal>(3, triangles);
	}

	finish("triangularize domain surface");
}

