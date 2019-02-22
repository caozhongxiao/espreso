
#include "meshpreprocessing.h"

#include "esinfo/ecfinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.h"

#include "basis/containers/point.h"
#include "basis/containers/serializededata.h"
#include "basis/matrices/denseMatrix.h"
#include "basis/utilities/communication.h"
#include "basis/utilities/utils.h"

#include "mesh/mesh.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementsregionstore.h"
#include "mesh/store/boundaryregionstore.h"
#include "mesh/store/surfacestore.h"
#include "mesh/elements/element.h"

#include <algorithm>
#include <numeric>

using namespace espreso;

MeshPreprocessing::MeshPreprocessing(Mesh *mesh)
: _mesh(mesh), _morphing(NULL)
{

}

MeshPreprocessing::~MeshPreprocessing()
{

}

void MeshPreprocessing::linkNodesAndElements()
{
	linkNodesAndElements(
			_mesh->nodes->elements,
			_mesh->elements->procNodes,
			_mesh->elements->IDs,
			_mesh->elements->distribution);
}

void MeshPreprocessing::linkNodesAndElements(
		serializededata<esint, esint>* &nelements,
		serializededata<esint, esint> *enodes,
		serializededata<esint, esint> *eIDs,
		std::vector<size_t> &edistribution)
{
	eslog::startln("MESH: LINK NODES AND ELEMENTS");

	size_t threads = info::env::OMP_NUM_THREADS;

	serializededata<esint, esint> *nIDs = _mesh->nodes->IDs;
	serializededata<esint, esint> *nranks = _mesh->nodes->ranks;
	std::vector<size_t> &ndistribution = _mesh->nodes->distribution;

	// thread x neighbor x vector(from, to)
	std::vector<std::vector<std::vector<std::pair<esint, esint> > > > sBuffer(threads);
	std::vector<std::vector<std::pair<esint, esint> > > rBuffer(_mesh->neighbours.size());
	std::vector<std::pair<esint, esint> > localLinks;

	localLinks.resize(enodes->cend()->begin() - enodes->cbegin()->begin());
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto tnodes = enodes->cbegin(t);
		size_t offset = enodes->cbegin(t)->begin() - enodes->cbegin()->begin();

		for (size_t e = edistribution[t]; e < edistribution[t + 1]; ++e, ++tnodes) {
			for (auto n = tnodes->begin(); n != tnodes->end(); ++n, ++offset) {
				localLinks[offset].first = *n;
				localLinks[offset].second = eIDs->datatarray()[e];
			}
		}
	}

	utils::sortWithInplaceMerge(localLinks, enodes->datatarray().distribution());

	std::vector<size_t> tbegin(threads);
	for (size_t t = 1; t < threads; t++) {
		tbegin[t] = std::lower_bound(localLinks.begin() + tbegin[t - 1], localLinks.end(), ndistribution[t], [] (std::pair<esint, esint> &p, esint n) { return p.first < n; }) - localLinks.begin();
	}

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto ranks = nranks->cbegin(t);
		std::vector<std::vector<std::pair<esint, esint> > > tBuffer(_mesh->neighbours.size());

		auto begin = localLinks.begin() + tbegin[t];
		auto end = begin;

		auto send = [&] (esint id) {
			for (auto it = begin; it != end; ++it) {
				it->first = id;
			}
			size_t i = 0;
			for (auto rank = ranks->begin(); rank != ranks->end(); ++rank) {
				if (*rank != info::mpi::rank) {
					while (_mesh->neighbours[i] < *rank) ++i;
					tBuffer[i].insert(tBuffer[i].end(), begin, end);
				}
			}
		};

		for (size_t n = ndistribution[t]; n < ndistribution[t + 1]; ++n, ++ranks) {
			while (begin != localLinks.end() && begin->first < (esint)n) ++begin;
			if (begin != localLinks.end() && begin->first == (esint)n) {
				end = begin;
				while (end != localLinks.end() && end->first == begin->first) ++end;
				send(nIDs->datatarray()[n]);
				begin = end;
			}
		}

		sBuffer[t].swap(tBuffer);
	}

	for (size_t t = 1; t < threads; t++) {
		for (size_t r = 0; r < sBuffer[0].size(); r++) {
			sBuffer[0][r].insert(sBuffer[0][r].end(), sBuffer[t][r].begin(), sBuffer[t][r].end());
		}
	}

	if (!Communication::exchangeUnknownSize(sBuffer[0], rBuffer, _mesh->neighbours)) {
		eslog::error("ESPRESO internal error: addLinkFromTo - exchangeUnknownSize.\n");
	}

	std::vector<size_t> boundaries = { 0, localLinks.size() };
	for (size_t r = 0; r < rBuffer.size(); r++) {
		localLinks.insert(localLinks.end(), rBuffer[r].begin(), rBuffer[r].end());
		boundaries.push_back(localLinks.size());
	}

	utils::mergeAppendedData(localLinks, boundaries);

	std::vector<std::vector<esint> > linksBoundaries(threads);
	std::vector<std::vector<esint> > linksData(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		if (ndistribution[t] != ndistribution[t + 1]) {
			auto llink = std::lower_bound(localLinks.begin(), localLinks.end(), nIDs->datatarray()[ndistribution[t]], [] (std::pair<esint, esint> &p, esint n) { return p.first < n; });
			esint current;

			std::vector<esint> tBoundaries, tData;
			if (t == 0) {
				tBoundaries.push_back(0);
			}

			for (size_t n = ndistribution[t]; n < ndistribution[t + 1] && llink != localLinks.end(); ++n) {
				current = llink->first;
				while (llink != localLinks.end() && current == llink->first) {
					tData.push_back(llink->second);
					++llink;
				}
				tBoundaries.push_back(llink - localLinks.begin());
			}

			linksBoundaries[t].swap(tBoundaries);
			linksData[t].swap(tData);
		}
	}

	nelements = new serializededata<esint, esint>(linksBoundaries, linksData);

	eslog::endln("MESH: NODES AND ELEMENTS LINKED");
}

void MeshPreprocessing::exchangeHalo()
{
	// halo elements are all elements that have some shared node
	if (_mesh->nodes->elements == NULL) {
		this->linkNodesAndElements();
	}

	if (_mesh->elements->regions == NULL) {
		fillRegionMask();
	}

	eslog::startln("MESH: EXCHANGE HALO");

	std::vector<esint> eDistribution = _mesh->elements->gatherElementsProcDistribution();

	size_t threads = info::env::OMP_NUM_THREADS;
	std::vector<std::vector<esint> > sBuffer(_mesh->neighbours.size()), rBuffer(_mesh->neighbours.size());

	std::vector<std::vector<std::vector<esint> > > hElements(threads);

	// we have to got through all nodes because intervals are not computed yet
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<std::vector<esint> > telements(_mesh->neighbours.size());
		auto elinks = _mesh->nodes->elements->cbegin(t);
		size_t i = 0;

		for (auto ranks = _mesh->nodes->ranks->cbegin(t); ranks != _mesh->nodes->ranks->cend(t); ++ranks, ++elinks) {
			auto begin = elinks->begin();
			auto end = elinks->begin();
			if (ranks->size() > 1) {
				i = 0;
				while (begin != elinks->end() && *begin < eDistribution[info::mpi::rank]) ++begin;
				end = begin;
				while (end != elinks->end() && *end < eDistribution[info::mpi::rank + 1]) ++end;
				for (auto rank = ranks->begin(); rank != ranks->end(); ++rank) {
					if (*rank != info::mpi::rank) {
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
		utils::sortWithInplaceMerge(hElements[0][n], tdist[n]);
	}
	#pragma omp parallel for
	for (size_t n = 0; n < _mesh->neighbours.size(); ++n) {
		utils::removeDuplicity(hElements[0][n]);
		tdist[n] = tarray<size_t>::distribute(threads, hElements[0][n].size());
		sBuffer[n].resize((4 + rsize) * hElements[0][n].size());
	}

	esint offset = eDistribution[info::mpi::rank];
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
				sBuffer[n][(4 + rsize) * e + 3] = (esint)code[hElements[0][n][e] - offset]->code;
				memcpy(sBuffer[n].data() + (4 + rsize) * e + 4, regions.data() + (hElements[0][n][e] - offset) * rsize, sizeof(esint) * rsize);
			}
		}
	}

	if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, _mesh->neighbours)) {
		eslog::error("ESPRESO internal error: exchange halo elements.\n");
	}

	std::vector<std::vector<esint> > hid(threads), hregions(threads);
	std::vector<std::vector<int> > hbody(threads), hmaterial(threads);

	std::vector<std::vector<Element*> > hcode(threads);

	for (size_t n = 0; n < rBuffer.size(); ++n) {
		std::vector<size_t> distribution = tarray<size_t>::distribute(threads, rBuffer[n].size() / (4 + rsize));
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (size_t e = distribution[t]; e < distribution[t + 1]; ++e) {
				hid[t].push_back(rBuffer[n][(4 + rsize) * e + 0]);
				hbody[t].push_back(rBuffer[n][(4 + rsize) * e + 1]);
				hmaterial[t].push_back(rBuffer[n][(4 + rsize) * e + 2]);
				hcode[t].push_back(_mesh->edata + rBuffer[n][(4 + rsize) * e + 3]);
				hregions[t].insert(hregions[t].end(), rBuffer[n].data() + (4 + rsize) * e + 4, rBuffer[n].data() + (4 + rsize) * e + 4 + rsize);
			}
		}
	}

	_mesh->halo->IDs = new serializededata<esint, esint>(1, hid);
	_mesh->halo->body = new serializededata<esint, int>(1, hbody);
	_mesh->halo->material = new serializededata<esint, int>(1, hmaterial);
	_mesh->halo->epointers = new serializededata<esint, Element*>(1, hcode);
	_mesh->halo->regions = new serializededata<esint, esint>(rsize, hregions);

	_mesh->halo->size = _mesh->halo->IDs->datatarray().size();
	_mesh->halo->distribution = _mesh->halo->IDs->datatarray().distribution();

	const auto &hIDs = _mesh->halo->IDs->datatarray();
	std::vector<esint> permutation(hIDs.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (esint i, esint j) { return hIDs[i] < hIDs[j]; });
	_mesh->halo->permute(permutation);

	eslog::endln("MESH: HALO EXCHANGED");
}

void MeshPreprocessing::computeElementsNeighbors()
{
	computeElementsNeighbors(
			_mesh->nodes->elements,
			_mesh->elements->neighbors,
			_mesh->elements->procNodes,
			_mesh->elements->IDs,
			_mesh->elements->epointers,
			_mesh->elements->distribution);
}

void MeshPreprocessing::computeElementsNeighbors(
		serializededata<esint, esint>* &nelements,
		serializededata<esint, esint>* &eneighbors,
		serializededata<esint, esint> *enodes,
		serializededata<esint, esint> *eIDs,
		serializededata<esint, Element*> *epointers,
		std::vector<size_t> &edistribution)
{
	if (nelements == NULL) {
		this->linkNodesAndElements(nelements, enodes, eIDs, edistribution);
	}

	eslog::startln("MESH: COMPUTE ELEMENTS NEIGHBOURS");

	size_t threads = info::env::OMP_NUM_THREADS;

	std::vector<std::vector<esint> > dualDistribution(threads);
	std::vector<std::vector<esint> > dualData(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto nodes = enodes->cbegin(t);

		std::vector<esint> tdist, tdata, intersection;
		if (t == 0) {
			tdist.push_back(0);
		}

		for (size_t e = edistribution[t]; e < edistribution[t + 1]; ++e, ++nodes) {
			for (auto face = epointers->datatarray()[e]->faces->begin(); face != epointers->datatarray()[e]->faces->end(); ++face) {
				auto telements = nelements->cbegin() + nodes->at(*face->begin());
				intersection.clear();
				for (auto n = telements->begin(); n != telements->end(); ++n) {
					if (*n != eIDs->datatarray()[e]) {
						intersection.push_back(*n);
					}
				}
				for (auto n = face->begin() + 1; n != face->end() && intersection.size(); ++n) {
					telements = nelements->cbegin() + nodes->at(*n);
					auto it1 = intersection.begin();
					auto it2 = telements->begin();
					auto last = intersection.begin();
					while (it1 != intersection.end()) {
						while (it2 != telements->end() && *it2 < *it1) {
							++it2;
						}
						if (it2 == telements->end()) {
							break;
						}
						if (*it1 == *it2) {
							*last++ = *it1++;
						} else {
							it1++;
						}
					}
					intersection.resize(last - intersection.begin());
				}
				tdata.push_back(intersection.size() ? intersection.front() : -1);
			}
			tdist.push_back(tdata.size());
		}

		dualDistribution[t].swap(tdist);
		dualData[t].swap(tdata);
	}

	utils::threadDistributionToFullDistribution(dualDistribution);

	eneighbors = new serializededata<esint, esint>(dualDistribution, dualData);

	eslog::endln("MESH: ELEMENTS NEIGHBOURS COMPUTED");
}

void MeshPreprocessing::computeSurfaceElementNeighbors(SurfaceStore *surface)
{
	surface->eoffset = surface->elements->structures();
	surface->IDs = new serializededata<esint, esint>(1, tarray<esint>(info::env::OMP_NUM_THREADS, surface->eoffset));
	Communication::exscan(surface->eoffset);
	std::iota(surface->IDs->datatarray().begin(), surface->IDs->datatarray().end(), surface->eoffset);
	computeElementsNeighbors(
			surface->nelements,
			surface->neighbors,
			surface->elements,
			surface->IDs,
			surface->epointers,
			surface->edistribution);
}

void MeshPreprocessing::computeElementsCenters()
{
	eslog::startln("MESH: COMPUTE ELEMENTS CENTERS");

	size_t threads = info::env::OMP_NUM_THREADS;

	_mesh->elements->centers = new serializededata<esint, float>(_mesh->dimension, tarray<float>(threads, _mesh->dimension * _mesh->elements->size));

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		Point center;
		size_t eindex = _mesh->elements->distribution[t];
		for (auto e = _mesh->elements->procNodes->cbegin(t); e != _mesh->elements->procNodes->cend(t); ++e, ++eindex) {
			center.x = center.y = center.z = 0;
			for (auto n = e->begin(); n != e->end(); ++n) {
				center += _mesh->nodes->coordinates->datatarray()[*n];
			}
			center /= e->size();
			_mesh->elements->centers->datatarray()[_mesh->dimension * eindex + 0] = center.x;
			_mesh->elements->centers->datatarray()[_mesh->dimension * eindex + 1] = center.y;
			if (_mesh->dimension == 3) {
				_mesh->elements->centers->datatarray()[_mesh->dimension * eindex + 2] = center.z;
			}
		}
	}
	eslog::endln("MESH: ELEMENTS CENTERS COMPUTED");
}

void MeshPreprocessing::computeDecomposedDual(std::vector<esint> &dualDist, std::vector<esint> &dualData)
{
	bool separateRegions = info::ecf->decomposition.separate_regions;
	bool separateMaterials = info::ecf->decomposition.separate_materials;
	bool separateEtypes = info::ecf->decomposition.separate_etypes;

	if (_mesh->elements->neighbors == NULL) {
		this->computeElementsNeighbors();
	}

	if (separateRegions && _mesh->elements->regions == NULL) {
		this->fillRegionMask();
	}

	eslog::startln("MESH: COMPUTE LOCAL DUAL GRAPH");

	size_t threads = info::env::OMP_NUM_THREADS;
	esint eBegin = _mesh->elements->gatherElementsProcDistribution()[info::mpi::rank];
	esint eEnd   = eBegin + _mesh->elements->size;

	std::vector<esint> dDistribution(_mesh->elements->size + 1);
	std::vector<std::vector<esint> > dData(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<esint> tdata;
		int mat1 = 0, mat2 = 0, reg = 0, etype1 = 0, etype2 = 0;
		int rsize = _mesh->elements->regionMaskSize;

		auto neighs = _mesh->elements->neighbors->cbegin(t);
		for (size_t e = _mesh->elements->distribution[t]; e < _mesh->elements->distribution[t + 1]; ++e, ++neighs) {
			for (auto n = neighs->begin(); n != neighs->end(); ++n) {
				if (*n != -1 && eBegin <= *n && *n < eEnd) {
					if (separateMaterials) {
						mat1 = _mesh->elements->material->datatarray()[e];
						mat2 = _mesh->elements->material->datatarray()[*n - eBegin];
					}
					if (separateRegions) {
						reg = memcmp(_mesh->elements->regions->datatarray().data() + e * rsize, _mesh->elements->regions->datatarray().data() + (*n - eBegin) * rsize, sizeof(esint) * rsize);
					}
					if (separateEtypes) {
						etype1 = (int)_mesh->elements->epointers->datatarray()[e]->type;
						etype2 = (int)_mesh->elements->epointers->datatarray()[*n - eBegin]->type;
					}

					if (mat1 == mat2 && !reg && etype1 == etype2) {
						tdata.push_back(*n - eBegin);
					}
				}
			}
			dDistribution[e + 1] = tdata.size();
		}

		dData[t].swap(tdata);
	}

	utils::threadDistributionToFullDistribution(dDistribution, _mesh->elements->distribution);
	for (size_t t = 1; t < threads; t++) {
		dData[0].insert(dData[0].end(), dData[t].begin(), dData[t].end());
	}

	dualDist.swap(dDistribution);
	dualData.swap(dData[0]);

	eslog::endln("MESH: LOCAL DUAL GRAPH COMPUTED");
}

void MeshPreprocessing::computeRegionsSurface()
{
	if (_mesh->elements->neighbors == NULL) {
		this->computeElementsNeighbors();
	}
	if (_mesh->halo->IDs == NULL) {
		this->exchangeHalo();
	}

	eslog::startln("MESH: COMPUTE REGION SURFACE");

	size_t threads = info::env::OMP_NUM_THREADS;
	esint eBegin = _mesh->elements->gatherElementsProcDistribution()[info::mpi::rank];
	esint eEnd = eBegin + _mesh->elements->size;

	for (size_t r = 0; r < _mesh->elementsRegions.size(); r++) {
		std::vector<std::vector<esint> > faces(threads), facesDistribution(threads), ecounters(threads, std::vector<esint>((int)Element::CODE::SIZE));
		std::vector<std::vector<Element*> > fpointers(threads);

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			esint hindex, addFace = 0;
			int rsize = _mesh->elements->regionMaskSize;
			auto nodes = _mesh->elements->procNodes->cbegin();
			auto neighs = _mesh->elements->neighbors->cbegin();
			const auto &regions = _mesh->elements->regions->datatarray();
			const auto &epointers = _mesh->elements->epointers->datatarray();

			std::vector<esint> fdist, fdata, ecounter((int)Element::CODE::SIZE);
			std::vector<Element*> fpointer;
			if (t == 0) {
				fdist.push_back(0);
			}

			esint prev = 0;
			for (auto e = _mesh->elementsRegions[r]->elements->datatarray().cbegin(t); e != _mesh->elementsRegions[r]->elements->datatarray().cend(t); prev = *e++) {
				nodes += *e - prev;
				neighs += *e - prev;
				for (size_t n = 0; n < neighs->size(); ++n) {
					if (neighs->at(n) != -1 && r) {
						if (neighs->at(n) < eBegin || eEnd <= neighs->at(n)) {
							hindex = std::lower_bound(_mesh->halo->IDs->datatarray().begin(), _mesh->halo->IDs->datatarray().end(), neighs->at(n)) - _mesh->halo->IDs->datatarray().begin();
							addFace = memcmp(regions.data() + *e * rsize, _mesh->halo->regions->datatarray().data() + hindex * rsize, sizeof(esint) * rsize);
						} else {
							addFace = memcmp(regions.data() + *e * rsize, regions.data() + (neighs->at(n) - eBegin) * rsize, sizeof(esint) * rsize);
						}
					} else {
						addFace = neighs->at(n) == -1;
					}
					if (addFace) {
						auto face = epointers[*e]->faces->begin() + n;
						for (auto f = face->begin(); f != face->end(); ++f) {
							fdata.push_back(nodes->at(*f));
						}
						fdist.push_back(fdata.size());
						fpointer.push_back(epointers[*e]->facepointers->datatarray()[n]);
						++ecounter[(int)fpointer.back()->code];
						addFace = 0;
					}
				}
			}

			facesDistribution[t].swap(fdist);
			faces[t].swap(fdata);
			fpointers[t].swap(fpointer);
			ecounters[t].swap(ecounter);
		}

		if (_mesh->elementsRegions[r]->surface == NULL) {
			_mesh->elementsRegions[r]->surface = new SurfaceStore();
		}

		for (size_t t = 1; t < threads; t++) {
			for (size_t e = 0; e < ecounters[0].size(); e++) {
				ecounters[0][e] += ecounters[t][e];
			}
		}

		serializededata<esint, Element*>::balance(1, fpointers);
		_mesh->elementsRegions[r]->surface->epointers = new serializededata<esint, Element*>(1, fpointers);
		_mesh->elementsRegions[r]->surface->ecounters = ecounters[0];

		_mesh->elementsRegions[r]->surface->edistribution = _mesh->elementsRegions[r]->surface->epointers->datatarray().distribution();

		if (
				_mesh->elementsRegions[r]->surface->edistribution.back() &&
				_mesh->elementsRegions[r]->surface->ecounters[(int)Element::CODE::TRIANGLE3] == (esint)_mesh->elementsRegions[r]->surface->edistribution.back()) {

			serializededata<esint, esint>::balance(3, faces, &_mesh->elementsRegions[r]->surface->edistribution);
			_mesh->elementsRegions[r]->surface->elements = new serializededata<esint, esint>(3, faces);
			_mesh->elementsRegions[r]->surface->triangles = _mesh->elementsRegions[r]->surface->elements;
			_mesh->elementsRegions[r]->surface->tdistribution = _mesh->elementsRegions[r]->surface->edistribution;
		} else {
			utils::threadDistributionToFullDistribution(facesDistribution);
			serializededata<esint, esint>::balance(facesDistribution, faces, &_mesh->elementsRegions[r]->surface->edistribution);
			_mesh->elementsRegions[r]->surface->elements = new serializededata<esint, esint>(facesDistribution, faces);
		}
	}

	eslog::endln("MESH: REGION SURFACE COMPUTED");
}

void MeshPreprocessing::triangularizeSurface(SurfaceStore *surface)
{
	if (surface == NULL) {
		return;
	}

	eslog::startln("MESH: TRIANGULARIZE SURFACE");

	size_t threads = info::env::OMP_NUM_THREADS;

	if (surface->triangles == NULL) {

		std::vector<std::vector<esint> > triangles(threads);
		std::vector<std::vector<size_t> > intervals(threads);


		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			std::vector<esint> ttriangles;
			std::vector<size_t> tintervals;
			if (t == 0) {
				tintervals.push_back(0);
			}

			auto elements = surface->elements->cbegin(t);
			const auto &epointers = surface->epointers->datatarray().begin();

			for (size_t e = surface->edistribution[t]; e < surface->edistribution[t + 1]; ++e, ++elements) {
				for (auto n = epointers[e]->triangles->datatarray().cbegin(); n != epointers[e]->triangles->datatarray().cend(); ++n) {
					ttriangles.push_back(elements->at(*n));
				}
			}
			tintervals.push_back(ttriangles.size() / 3);

			intervals[t].swap(tintervals);
			triangles[t].swap(ttriangles);
		}

		utils::threadDistributionToFullDistribution(intervals);
		utils::mergeThreadedUniqueData(intervals);

		surface->tdistribution = intervals[0];
		surface->triangles = new serializededata<esint, esint>(3, triangles);
	}

	eslog::endln("MESH: SURFACE TRIANGULARIZED");
}

void MeshPreprocessing::triangularizeBoundary(BoundaryRegionStore *boundary)
{
	eslog::startln("MESH: TRIANGULARIZE BOUNDARY");

	size_t threads = info::env::OMP_NUM_THREADS;

	if (boundary->dimension == 2) {

		std::vector<std::vector<esint> > triangles(threads);
		std::vector<std::vector<size_t> > intervals(threads);


		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			std::vector<esint> ttriangles;
			std::vector<size_t> tintervals;
			if (t == 0) {
				tintervals.push_back(0);
			}

			auto elements = boundary->procNodes->cbegin(t);
			const auto &epointers = boundary->epointers->datatarray().begin();

			for (size_t e = boundary->distribution[t]; e < boundary->distribution[t + 1]; ++e, ++elements) {
				for (auto n = epointers[e]->triangles->datatarray().cbegin(); n != epointers[e]->triangles->datatarray().cend(); ++n) {
					ttriangles.push_back(elements->at(*n));
				}
			}
			tintervals.push_back(ttriangles.size() / 3);

			intervals[t].swap(tintervals);
			triangles[t].swap(ttriangles);
		}

		utils::threadDistributionToFullDistribution(intervals);
		utils::mergeThreadedUniqueData(intervals);

		boundary->triangles = new serializededata<esint, esint>(3, triangles);
	}

	eslog::endln("MESH: BOUNDARY TRIANGULARIZED");
}

void MeshPreprocessing::computeBoundaryNodes(std::vector<esint> &externalBoundary, std::vector<esint> &internalBoundary)
{
	if (_mesh->elements->neighbors == NULL) {
		this->computeElementsNeighbors();
	}

	eslog::startln("MESH: COMPUTE BOUNDARY NODES");

	size_t threads = info::env::OMP_NUM_THREADS;

	std::vector<std::vector<esint> > external(threads), internal(threads);

	esint eoffset = _mesh->elements->gatherElementsProcDistribution()[info::mpi::rank];

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<esint> texternal, tinternal;

		auto neighbors = _mesh->elements->neighbors->cbegin() + _mesh->elements->elementsDistribution[_mesh->elements->domainDistribution[t]];
		auto enodes = _mesh->elements->procNodes->cbegin() + _mesh->elements->elementsDistribution[_mesh->elements->domainDistribution[t]];
		for (esint d = _mesh->elements->domainDistribution[t]; d < _mesh->elements->domainDistribution[t + 1]; d++) {
			esint dbegin = _mesh->elements->elementsDistribution[d];
			esint dend = _mesh->elements->elementsDistribution[d + 1];

			for (esint e = dbegin; e < dend; ++e, ++neighbors, ++enodes) {
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

	utils::sortWithUniqueMerge(internal);
	utils::sortWithUniqueMerge(external);

	externalBoundary.swap(external[0]);

	auto n2i = [ & ] (size_t neighbour) {
		return std::lower_bound(_mesh->neighbours.begin(), _mesh->neighbours.end(), neighbour) - _mesh->neighbours.begin();
	};

	// external nodes need to be synchronized
	std::vector<std::vector<esint> > sBuffer(_mesh->neighbours.size()), rBuffer(_mesh->neighbours.size());
	std::vector<esint> nExternal;

	for (size_t i = 0; i < externalBoundary.size(); i++) {
		auto nrank = _mesh->nodes->ranks->cbegin() + externalBoundary[i];
		for (auto rank = nrank->begin(); rank != nrank->end(); ++rank) {
			if (*rank != info::mpi::rank) {
				sBuffer[n2i(*rank)].push_back(_mesh->nodes->IDs->datatarray()[externalBoundary[i]]);
			}
		}
	}

	for (size_t n = 0; n < _mesh->neighbours.size(); n++) {
		std::sort(sBuffer[n].begin(), sBuffer[n].end());
	}

	if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, _mesh->neighbours)) {
		eslog::error("ESPRESO internal error: exchange external nodes.\n");
	}

	for (size_t n = 0; n < _mesh->neighbours.size(); n++) {
		nExternal.insert(nExternal.end(), rBuffer[n].begin(), rBuffer[n].end());
	}
	utils::sortAndRemoveDuplicity(nExternal);

	for (size_t n = 0; n < _mesh->neighbours.size(); n++) {
		nExternal.resize(std::set_difference(nExternal.begin(), nExternal.end(), sBuffer[n].begin(), sBuffer[n].end(), nExternal.begin()) - nExternal.begin());
	}

	for (size_t n = 0; n < nExternal.size(); n++) {
		std::vector<std::vector<esint> > tnExternal(threads);
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

	eslog::endln("MESH: BOUNDARY NODES COMPUTED");
}

void MeshPreprocessing::computeRegionArea(BoundaryRegionStore *store)
{
	double A = 0;
	auto nodes = store->procNodes->cbegin();
	const auto &epointers = store->epointers->datatarray();
	const auto &coordinates = _mesh->nodes->coordinates->datatarray();
	for (size_t e = 0; e < store->procNodes->structures(); ++e, ++nodes) {

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

	MPI_Allreduce(&A, &store->area, 1, MPI_DOUBLE, MPI_SUM, info::mpi::comm);
}
