
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
#include "../../basis/utilities/communication.h"
#include "../../basis/utilities/utils.h"
#include "../../basis/utilities/parser.h"
#include "../../basis/logging/logging.h"
#include "../../basis/logging/timeeval.h"

#include "../../config/ecf/root.h"
#include "../../config/ecf/decomposition.h"

#include <algorithm>
#include <numeric>
#include <cstring>
#include "../../wrappers/metis/wmetis.h"
#include "../../wrappers/metis/wparmetis.h"

using namespace espreso;


void MeshPreprocessing::reclusterize()
{
	if (environment->MPIsize == 1) {
		skip("re-distribution of the mesh to processes");
		return;
	}

	start("re-distribution of the mesh to processes");

	if (_mesh->elements->neighbors == NULL) {
		this->computeElementsNeighbors();
	}

	if (_mesh->halo->IDs == NULL) {
		exchangeHalo();
	}

	bool separateRegions = _mesh->configuration.decomposition.separate_regions;
	bool separateMaterials = _mesh->configuration.decomposition.separate_materials;
	bool separateEtypes = _mesh->configuration.decomposition.separate_etypes;

	if (separateRegions && _mesh->elements->regions == NULL) {
		fillRegionMask();
	}

	size_t threads = environment->OMP_NUM_THREADS;

	eslocal eoffset = _mesh->elements->gatherElementsProcDistribution()[environment->MPIrank];

	std::vector<eslocal> dDistribution(_mesh->elements->size + 1);
	std::vector<std::vector<eslocal> > dData(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<eslocal> tdata;
		int mat1 = 0, mat2 = 0, reg = 0, etype1 = 0, etype2 = 0;
		int rsize = _mesh->elements->regionMaskSize;
		eslocal hindex;

		auto neighs = _mesh->elements->neighbors->cbegin(t);
		for (size_t e = _mesh->elements->distribution[t]; e < _mesh->elements->distribution[t + 1]; ++e, ++neighs) {
			for (auto n = neighs->begin(); n != neighs->end(); ++n) {
				if (*n != -1) {
					if ((separateRegions || separateMaterials || separateEtypes) && (*n < eoffset || eoffset + _mesh->elements->size <= *n)) {
						hindex = std::lower_bound(_mesh->halo->IDs->datatarray().begin(), _mesh->halo->IDs->datatarray().end(), *n) - _mesh->halo->IDs->datatarray().begin();
					}
					if (separateMaterials) {
						mat1 = _mesh->elements->material->datatarray()[e];
						if (*n < eoffset || eoffset + _mesh->elements->size <= *n) {
							mat2 = _mesh->halo->material->datatarray()[hindex];
						} else {
							mat2 = _mesh->elements->material->datatarray()[*n - eoffset];
						}
					}
					if (separateRegions) {
						if (*n < eoffset || eoffset + _mesh->elements->size <= *n) {
							reg = memcmp(_mesh->elements->regions->datatarray().data() + e * rsize, _mesh->halo->regions->datatarray().data() + hindex * rsize, sizeof(int) * rsize);
						} else {
							reg = memcmp(_mesh->elements->regions->datatarray().data() + e * rsize, _mesh->elements->regions->datatarray().data() + (*n - eoffset) * rsize, sizeof(int) * rsize);
						}
					}
					if (separateEtypes) {
						etype1 = (int)_mesh->elements->epointers->datatarray()[e]->type;
						if (*n < eoffset || eoffset + _mesh->elements->size <= *n) {
							etype2 = (int)_mesh->halo->epointers->datatarray()[hindex]->type;
						} else {
							etype2 = (int)_mesh->elements->epointers->datatarray()[*n - eoffset]->type;
						}
					}

					if (mat1 == mat2 && !reg && etype1 == etype2) {
						tdata.push_back(*n);
					}
				}
			}
			dDistribution[e + 1] = tdata.size();
		}

		dData[t].swap(tdata);
	}

	Esutils::threadDistributionToFullDistribution(dDistribution, _mesh->elements->distribution);
	for (size_t t = 1; t < threads; t++) {
		dData[0].insert(dData[0].end(), dData[t].begin(), dData[t].end());
	}

	std::vector<eslocal> edistribution = _mesh->elements->gatherElementsProcDistribution();
	std::vector<eslocal> partition(_mesh->elements->size), permutation(_mesh->elements->size);
	start("ParMETIS::KWay");
	eslocal edgecut = ParMETIS::call(ParMETIS::METHOD::ParMETIS_V3_PartKway,
		edistribution.data(),
		dDistribution.data(), dData[0].data(),
		0, NULL,
		0, NULL, NULL,
		partition.data()
	);
	finish("ParMETIS::KWay");

	if (_mesh->configuration.decomposition.metis_options.adaptive_refinement) {
		start("ParMETIS::AdaptiveRepart");
		eslocal prev = 2 * edgecut;
		while (1.01 * edgecut < prev) {
			prev = edgecut;
			edgecut = ParMETIS::call(ParMETIS::METHOD::ParMETIS_V3_AdaptiveRepart,
				edistribution.data(),
				dDistribution.data(), dData[0].data(),
				0, NULL, // 3, _mesh->elements->coordinates->datatarray(),
				0, NULL, NULL,
				partition.data()
			);
		}
		finish("ParMETIS::AdaptiveRepart");
	}

	this->exchangeElements(partition);

	finish("re-distribution of the mesh to processes");
}

void MeshPreprocessing::partitiate(eslocal parts)
{
	start("decomposition of the mesh");

	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<eslocal> dualDist, dualData;
	this->computeDecomposedDual(dualDist, dualData);

	std::vector<int> partID(_mesh->elements->size, -1);

	int nextID = 0;
	for (eslocal e = 0; e < _mesh->elements->size; ++e) {
		std::vector<eslocal> stack;
		if (partID[e] == -1) {
			stack.push_back(e);
			partID[e] = nextID;
			while (stack.size()) {
				eslocal current = stack.back();
				stack.pop_back();
				for (auto n = dualData.begin() + dualDist[current]; n != dualData.begin() + dualDist[current + 1]; ++n) {
					if (partID[*n] == -1) {
						stack.push_back(*n);
						partID[*n] = nextID;
					}
				}
			}
			nextID++;
		}
	}

	std::vector<int> clusters;
	std::vector<eslocal> partition(_mesh->elements->size);
	if (nextID == 1) {
		start("METIS::KWay");
		METIS::call(
				_mesh->configuration.decomposition.metis_options,
				_mesh->elements->size,
				dualDist.data(), dualData.data(),
				0, NULL, NULL,
				parts, partition.data());
		finish("METIS::KWay");
		clusters.resize(parts, 0);
		_mesh->elements->nclusters = 1;

	} else { // non-continuous dual graph
		// thread x part x elements
		std::vector<std::vector<std::vector<eslocal> > > tdecomposition(threads, std::vector<std::vector<eslocal> >(nextID));
		std::vector<std::vector<eslocal> > tdualsize(threads, std::vector<eslocal>(nextID));

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (size_t e = _mesh->elements->distribution[t]; e < _mesh->elements->distribution[t + 1]; ++e) {
				tdecomposition[t][partID[e]].push_back(e);
				tdualsize[t][partID[e]] += dualDist[e + 1] - dualDist[e];
			}
		}
		std::vector<std::vector<eslocal> > foffsets(nextID), noffsets(nextID);
		std::vector<eslocal> partoffset(nextID);
		#pragma omp parallel for
		for (int p = 0; p < nextID; p++) {
			foffsets[p].push_back(0);
			noffsets[p].push_back(0);
			for (size_t t = 1; t < threads; t++) {
				foffsets[p].push_back(tdecomposition[0][p].size());
				noffsets[p].push_back(tdualsize[0][p]);
				tdecomposition[0][p].insert(tdecomposition[0][p].end(), tdecomposition[t][p].begin(), tdecomposition[t][p].end());
				tdualsize[0][p] += tdualsize[t][p];
			}
		}
		for (int p = 1; p < nextID; p++) {
			partoffset[p] = partoffset[p - 1] + tdecomposition[0][p - 1].size();
		}

		std::vector<std::vector<eslocal> > frames(nextID), neighbors(nextID);
		#pragma omp parallel for
		for (int p = 0; p < nextID; p++) {
			frames[p].resize(1 + tdecomposition[0][p].size());
			neighbors[p].resize(tdualsize[0][p]);
		}

		// TODO: try parallelization
		// #pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			size_t partindex;
			std::vector<eslocal> foffset(nextID), noffset(nextID);
			for (int p = 0; p < nextID; p++) {
				foffset[p] = foffsets[p][t];
				noffset[p] = noffsets[p][t];
			}

			for (size_t e = _mesh->elements->distribution[t]; e < _mesh->elements->distribution[t + 1]; ++e) {
				partindex = partID[e];

				frames[partindex][++foffset[partindex]] = dualDist[e + 1] - dualDist[e];
				if (e > _mesh->elements->distribution[t]) {
					frames[partindex][foffset[partindex]] += frames[partindex][foffset[partindex] - 1];
				} else {
					frames[partindex][foffset[partindex]] += noffset[partindex];
				}
				auto node = dualData.begin() + dualDist[e];
				for (eslocal n = frames[partindex][foffset[partindex]] - (dualDist[e + 1] - dualDist[e]); n < frames[partindex][foffset[partindex]]; ++n, ++node) {
					neighbors[partindex][n] = std::lower_bound(tdecomposition[0][partindex].begin(), tdecomposition[0][partindex].end(), *node) - tdecomposition[0][partindex].begin();
				}
			}
		}

		std::vector<eslocal> pparts(nextID);

		double averageDomainSize = _mesh->elements->size / (double)parts;
		size_t partsCounter = 0;
		for (int p = 0; p < nextID; p++) {
			partsCounter += pparts[p] = std::ceil((frames[p].size() - 1) / averageDomainSize);
			clusters.resize(partsCounter, p);
		}
		_mesh->elements->nclusters = nextID;

		start("METIS::KWay");
		#pragma omp parallel for
		for (int p = 0; p < nextID; p++) {
			METIS::call(
					_mesh->configuration.decomposition.metis_options,
					frames[p].size() - 1,
					frames[p].data(), neighbors[p].data(),
					0, NULL, NULL,
					pparts[p], partition.data() + partoffset[p]);
		}
		finish("METIS::KWay");

		std::vector<eslocal> ppartition = partition;
		nextID = 0;
		for (size_t p = 0; p < tdecomposition[0].size(); p++) {
			for (size_t i = 0; i < tdecomposition[0][p].size(); ++i) {
				partition[tdecomposition[0][p][i]] = ppartition[partoffset[p] + i] + nextID;
			}
			nextID += pparts[p];
		}
	}

	std::vector<eslocal> permutation(partition.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) {
		if (partition[i] == partition[j]) {
			return i < j;
		}
		return partition[i] < partition[j];
	});

	std::vector<size_t> domainDistribution;
	std::vector<size_t> tdistribution;

	eslocal partindex = 0;
	auto begin = permutation.begin();
	while (begin != permutation.end()) {
		domainDistribution.push_back(begin - permutation.begin());
		begin = std::lower_bound(begin, permutation.end(), ++partindex, [&] (eslocal i, eslocal val) {
			return partition[i] < val;
		});
	}
	domainDistribution.push_back(permutation.size());

	// TODO: improve domain distribution for more complicated decomposition
	if (domainDistribution.size() == threads + 1) {
		tdistribution = std::vector<size_t>(domainDistribution.begin(), domainDistribution.end());
	} else {
		if (domainDistribution.size() < threads + 1) {
			tdistribution = tarray<eslocal>::distribute(threads, permutation.size());
		} else {
			double averageThreadSize = _mesh->elements->size / (double)threads;
			tdistribution.push_back(0);
			for (size_t t = 0; t < threads - 1; t++) {
				auto more = std::lower_bound(domainDistribution.begin(), domainDistribution.end(), tdistribution.back() + averageThreadSize);
				auto less = more - 1;
				if (std::fabs(*less - averageThreadSize * (t + 1)) < std::fabs(*more - averageThreadSize * (t + 1))) {
					tdistribution.push_back(*less);
				} else {
					tdistribution.push_back(*more);
				}
			}
			tdistribution.push_back(permutation.size());
		}
	}

	for (size_t i = 1; i < domainDistribution.size(); i++) {
		if (domainDistribution[i - 1] != domainDistribution[i]) {
			_mesh->elements->clusters.push_back(clusters[i - 1]);
		}
	}
	Esutils::removeDuplicity(domainDistribution);

	std::vector<eslocal> domainCounter(threads);
	for (size_t t = 0; t < threads; t++) {
		if (domainDistribution.size() < threads + 1) {
			if (t < domainDistribution.size() - 1) {
				++domainCounter[t];
			}
		} else {
			auto begin = std::lower_bound(domainDistribution.begin(), domainDistribution.end(), tdistribution[t]);
			auto end   = std::lower_bound(domainDistribution.begin(), domainDistribution.end(), tdistribution[t + 1]);
			for (auto it = begin; it != end; ++it) {
				++domainCounter[t];
			}
		}
	}

	_mesh->elements->ndomains = Esutils::sizesToOffsets(domainCounter);
	_mesh->elements->firstDomain = _mesh->elements->ndomains;
	Communication::exscan(_mesh->elements->firstDomain, MPITools::operations().sizeToOffsetsEslocal);
	domainCounter.push_back(_mesh->elements->ndomains);
	_mesh->elements->domainDistribution = domainCounter;

	_mesh->elements->elementsDistribution.push_back(0);
	for (size_t t = 0; t < threads; t++) {
		if (domainDistribution.size() < threads + 1) {
			if (t < domainDistribution.size() - 1) {
				_mesh->elements->elementsDistribution.push_back(domainDistribution[t + 1]);
			}
		} else {
			auto begin = std::lower_bound(domainDistribution.begin(), domainDistribution.end(), tdistribution[t]);
			auto end   = std::lower_bound(domainDistribution.begin(), domainDistribution.end(), tdistribution[t + 1]);
			for (auto it = begin; it != end; ++it) {
				_mesh->elements->elementsDistribution.push_back(*(it + 1));
			}
		}
	}

	arrangeElementsPermutation(permutation);
	this->permuteElements(permutation, tdistribution);

	arrangeNodes();
//	arrangeRegions();

	finish("decomposition of the mesh");
}

void MeshPreprocessing::exchangeElements(const std::vector<eslocal> &partition)
{
	start("exchanging elements");

	// 0. Compute targets
	// 1. Serialize element data
	// 2. Serialize node data
	// 3. Send serialized data to target (new MPI processes)
	// 4. Deserialize data
	// 5. Balance nodes data to threads
	// 6. Re-index elements (IDs have to be always increasing)

	if (_mesh->nodes->elements == NULL) {
		// need for correctly update nodes ranks
		this->linkNodesAndElements();
	}

	size_t threads = environment->OMP_NUM_THREADS;

	TimeEval eval("exchange");
	eval.totalTime.startWithBarrier();

	TimeEvent e0("compute targets");
	e0.start();

	// Step 0: Compute targets

	std::vector<int> targets;
	{ // get target processes
		std::vector<std::vector<eslocal> > ttargets(threads);
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {

			std::vector<eslocal> ttt;
			std::vector<bool> tflags(environment->MPIsize, false);
			for (size_t e = _mesh->elements->distribution[t]; e < _mesh->elements->distribution[t + 1]; ++e) {
				if (partition[e] != environment->MPIrank) {
					tflags[partition[e]] = true;
				}
			}
			for (int r = 0; r < environment->MPIsize; r++) {
				if (tflags[r]) {
					ttt.push_back(r);
				}
			}

			ttargets[t].swap(ttt);
		}

		Esutils::mergeThreadedUniqueData(ttargets);
		targets = ttargets[0];
	}

	auto t2i = [ & ] (size_t target) {
		return std::lower_bound(targets.begin(), targets.end(), target) - targets.begin();
	};

	ElementStore *elements = new ElementStore(_mesh->_eclasses);

	std::vector<std::vector<eslocal> >  elemsIDs(threads);
	std::vector<std::vector<int> >      elemsBody(threads);
	std::vector<std::vector<int> >      elemsMaterial(threads);
	std::vector<std::vector<Element*> > elemsEpointer(threads);
	std::vector<std::vector<eslocal> >  elemsNodesDistribution(threads);
	std::vector<std::vector<eslocal> >  elemsNodesData(threads);
	std::vector<std::vector<eslocal> >  elemsRegions(threads);

	NodeStore *nodes = new NodeStore();

	std::vector<std::vector<eslocal> >  nodesIDs(threads);
	std::vector<std::vector<Point> >    nodesCoordinates(threads);
	std::vector<std::vector<eslocal> >  nodesElemsDistribution(threads);
	std::vector<std::vector<eslocal> >  nodesElemsData(threads);
	std::vector<std::vector<eslocal> >  nodesRegions(threads);

	std::vector<std::vector<std::vector<eslocal> > >  boundaryEDistribution(_mesh->boundaryRegions.size(), std::vector<std::vector<eslocal> >(threads));
	std::vector<std::vector<std::vector<eslocal> > >  boundaryEData(_mesh->boundaryRegions.size(), std::vector<std::vector<eslocal> >(threads));
	std::vector<std::vector<std::vector<Element*> > > boundaryEPointers(_mesh->boundaryRegions.size(), std::vector<std::vector<Element*> >(threads));

	// regions are transfered via mask
	int eregionsBitMaskSize = _mesh->elementsRegions.size() / (8 * sizeof(eslocal)) + (_mesh->elementsRegions.size() % (8 * sizeof(eslocal)) ? 1 : 0);
	int bregionsBitMaskSize = _mesh->boundaryRegions.size() / (8 * sizeof(eslocal)) + (_mesh->boundaryRegions.size() % (8 * sizeof(eslocal)) ? 1 : 0);

	// serialize data that have to be exchanged
	// the first thread value denotes the thread data size

	// threads x target x elements(id, body, material, code, dualsize, dualdata, nodesize, nodeindices)
	std::vector<std::vector<std::vector<eslocal> > > sElements(threads, std::vector<std::vector<eslocal> >(targets.size(), std::vector<eslocal>({ 0 })));
	std::vector<std::vector<eslocal> > rElements;

	// threads x target x nodes(id, point, linksize, links, regionMask) + size
	std::vector<std::vector<std::vector<eslocal> > > sNodes(threads, std::vector<std::vector<eslocal> >(targets.size(), std::vector<eslocal>({ 0 })));
	std::vector<std::vector<eslocal> > rNodes;

	// threads x target x boundary(prefix, (code, nodes))
	std::vector<std::vector<std::vector<eslocal> > > sBoundary(threads, std::vector<std::vector<eslocal> >(targets.size(), std::vector<eslocal>(_mesh->boundaryRegions.size())));
	std::vector<std::vector<eslocal> > rBoundary;

	e0.end();
	eval.addEvent(e0);

	TimeEvent e1("serialize elements");
	e1.start();

	// Step 1: Serialize element data

	std::vector<eslocal> regionElementMask(_mesh->elements->size * eregionsBitMaskSize);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		int maskOffset = 0;
		for (size_t r = 0; r < _mesh->elementsRegions.size(); r++) {
			maskOffset = r / (8 * sizeof(eslocal));
			int bit = 1 << (r % (8 * sizeof(eslocal)));
			auto begin = std::lower_bound(_mesh->elementsRegions[r]->elements->datatarray().begin(), _mesh->elementsRegions[r]->elements->datatarray().end(), _mesh->elements->distribution[t]);
			auto end = std::lower_bound(_mesh->elementsRegions[r]->elements->datatarray().begin(), _mesh->elementsRegions[r]->elements->datatarray().end(), _mesh->elements->distribution[t + 1]);
			for (auto i = begin; i != end; ++i) {
				regionElementMask[*i * eregionsBitMaskSize + maskOffset] |= bit;
			}
		}
	}

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto IDs = _mesh->elements->IDs->datatarray().data();
		auto body = _mesh->elements->body->datatarray().data();
		auto material = _mesh->elements->material->datatarray().data();
		auto code = _mesh->elements->epointers->datatarray().data();
		auto enodes = _mesh->elements->nodes->cbegin(t);
		auto nIDs = _mesh->nodes->IDs->datatarray().data();

		std::vector<std::vector<eslocal> > tsElements(targets.size(), std::vector<eslocal>({ 0 }));

		std::vector<eslocal>  telemsIDs;
		std::vector<int>      telemsBody;
		std::vector<int>      telemsMaterial;
		std::vector<Element*> telemsEpointer;
		std::vector<eslocal>  telemsNodesDistribution;
		std::vector<eslocal>  telemsNodesData;
		std::vector<eslocal>  telemsRegions;
		if (t == 0) {
			telemsNodesDistribution.push_back(0);
		}

		// estimation
		telemsIDs.reserve(1.5 * _mesh->elements->size / threads);
		telemsBody.reserve(1.5 * _mesh->elements->size / threads);
		telemsMaterial.reserve(1.5 * _mesh->elements->size / threads);
		telemsEpointer.reserve(1.5 * _mesh->elements->size / threads);
		telemsNodesDistribution.reserve(1.5 * _mesh->elements->size / threads);
		telemsRegions.reserve(1.5 * _mesh->elements->size / threads);

		size_t target;
		for (size_t e = _mesh->elements->distribution[t]; e < _mesh->elements->distribution[t + 1]; ++e, ++enodes) {
			if (partition[e] == environment->MPIrank) {
				telemsIDs.push_back(IDs[e]);
				telemsBody.push_back(body[e]);
				telemsMaterial.push_back(material[e]);
				telemsEpointer.push_back(code[e]);
				for (auto n = enodes->begin(); n != enodes->end(); ++n) {
					telemsNodesData.push_back(nIDs[*n]);
				}
				telemsNodesDistribution.push_back(telemsNodesData.size());
				telemsRegions.insert(telemsRegions.end(), regionElementMask.begin() + e * eregionsBitMaskSize, regionElementMask.begin() + (e + 1) * eregionsBitMaskSize);
			} else {
				target = t2i(partition[e]);
				tsElements[target].insert(tsElements[target].end(), { IDs[e], body[e], material[e], static_cast<eslocal>(code[e] - _mesh->_eclasses[t]) });
				tsElements[target].push_back(enodes->size());
				for (auto n = enodes->begin(); n != enodes->end(); ++n) {
					tsElements[target].push_back(nIDs[*n]);
				}
				tsElements[target].insert(tsElements[target].end(), regionElementMask.begin() + e * eregionsBitMaskSize, regionElementMask.begin() + (e + 1) * eregionsBitMaskSize);
			}
		}

		elemsIDs[t].swap(telemsIDs);
		elemsBody[t].swap(telemsBody);
		elemsMaterial[t].swap(telemsMaterial);
		elemsEpointer[t].swap(telemsEpointer);
		elemsNodesDistribution[t].swap(telemsNodesDistribution);
		elemsNodesData[t].swap(telemsNodesData);
		elemsRegions[t].swap(telemsRegions);

		sElements[t].swap(tsElements);
	}

	e1.end();
	eval.addEvent(e1);

	TimeEvent e2("serialize nodes");
	e2.start();

	// Step 2: Serialize node data

	std::vector<eslocal> regionNodeMask(_mesh->nodes->size * bregionsBitMaskSize);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		int maskOffset = 0;
		for (size_t r = 0; r < _mesh->boundaryRegions.size(); r++) {
			if (_mesh->boundaryRegions[r]->nodes) {
				maskOffset = r / (8 * sizeof(eslocal));
				int bit = 1 << (r % (8 * sizeof(eslocal)));
				auto begin = std::lower_bound(_mesh->boundaryRegions[r]->nodes->datatarray().begin(), _mesh->boundaryRegions[r]->nodes->datatarray().end(), _mesh->nodes->distribution[t]);
				auto end = std::lower_bound(_mesh->boundaryRegions[r]->nodes->datatarray().begin(), _mesh->boundaryRegions[r]->nodes->datatarray().end(), _mesh->nodes->distribution[t + 1]);
				for (auto i = begin; i != end; ++i) {
					regionNodeMask[*i * bregionsBitMaskSize + maskOffset] |= bit;
				}
			}
		}
	}

	eslocal eBegin = _mesh->elements->gatherElementsProcDistribution()[environment->MPIrank];
	eslocal eEnd = eBegin + _mesh->elements->size;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		const auto &IDs = _mesh->nodes->IDs->datatarray();
		const auto &coordinates = _mesh->nodes->coordinates->datatarray();
		auto elems = _mesh->nodes->elements->cbegin(t);

		std::vector<eslocal>  tnodesIDs;
		std::vector<Point>    tnodesCoordinates;
		std::vector<eslocal>  tnodesElemsDistribution;
		std::vector<eslocal>  tnodesElemsData;
		std::vector<eslocal>  tnodesRegions;

		if (t == 0) {
			tnodesElemsDistribution.push_back(0);
		}

		std::vector<std::vector<eslocal> > tsNodes(targets.size(), std::vector<eslocal>({ 0 }));

		tnodesIDs.reserve(1.5 * _mesh->nodes->size / threads);
		tnodesCoordinates.reserve(1.5 * _mesh->nodes->size / threads);
		tnodesElemsDistribution.reserve(1.5 * _mesh->nodes->size / threads);
		tnodesRegions.reserve(1.5 * _mesh->nodes->size / threads);

		size_t target;
		std::vector<bool> last(targets.size() + 1); // targets + me
		for (size_t n = _mesh->nodes->distribution[t]; n < _mesh->nodes->distribution[t + 1]; ++n, ++elems) {
			std::fill(last.begin(), last.end(), false);
			for (auto e = elems->begin(); e != elems->end(); ++e) {
				if (eBegin <= *e && *e < eEnd) {
					target = t2i(partition[*e - eBegin]);
					if (!last[target] && partition[*e - eBegin] != environment->MPIrank) {
						tsNodes[target].push_back(IDs[n]);
						tsNodes[target].insert(tsNodes[target].end(), reinterpret_cast<const eslocal*>(coordinates.data() + n), reinterpret_cast<const eslocal*>(coordinates.data() + n + 1));
						tsNodes[target].push_back(elems->size());
						tsNodes[target].insert(tsNodes[target].end(), elems->begin(), elems->end());
						tsNodes[target].insert(tsNodes[target].end(), regionNodeMask.begin() + n * bregionsBitMaskSize, regionNodeMask.begin() + (n + 1) * bregionsBitMaskSize);
						last[target] = true;
					}
					if (!last.back() && partition[*e - eBegin] == environment->MPIrank) {
						tnodesIDs.push_back(IDs[n]);
						tnodesCoordinates.push_back(coordinates[n]);
						tnodesElemsData.insert(tnodesElemsData.end(), elems->begin(), elems->end());
						tnodesElemsDistribution.push_back(tnodesElemsData.size());
						tnodesRegions.insert(tnodesRegions.end(), regionNodeMask.begin() + n * bregionsBitMaskSize, regionNodeMask.begin() + (n + 1) * bregionsBitMaskSize);
						last.back() = true;
					}
				}
			}
		}

		nodesIDs[t].swap(tnodesIDs);
		nodesCoordinates[t].swap(tnodesCoordinates);
		nodesElemsDistribution[t].swap(tnodesElemsDistribution);
		nodesElemsData[t].swap(tnodesElemsData);
		nodesRegions[t].swap(tnodesRegions);

		sNodes[t].swap(tsNodes);
	}

	e2.end();
	eval.addEvent(e2);

	TimeEvent e21("serialize boundary");
	e21.start();

	// Step 2.1: Serialize boundary regions data

	std::vector<eslocal> emembership;

	for (size_t r = 0; r < _mesh->boundaryRegions.size(); r++) {
		if (_mesh->boundaryRegions[r]->dimension) {
			emembership.clear();
			emembership.resize(_mesh->boundaryRegions[r]->distribution.back());
			std::vector<size_t> distribution = _mesh->boundaryRegions[r]->distribution;

			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				auto enodes = _mesh->boundaryRegions[r]->elements->cbegin() + distribution[t];
				std::vector<eslocal> nlinks;
				int counter;
				for (size_t e = distribution[t]; e < distribution[t + 1]; ++e, ++enodes) {
					nlinks.clear();
					for (auto n = enodes->begin(); n != enodes->end(); ++n) {
						auto links = _mesh->nodes->elements->cbegin() + *n;
						nlinks.insert(nlinks.end(), links->begin(), links->end());
					}
					std::sort(nlinks.begin(), nlinks.end());
					counter = 1;
					for (size_t i = 1; i < nlinks.size(); ++i) {
						if (nlinks[i - 1] == nlinks[i]) {
							++counter;
							if (counter == enodes->size() && eBegin <= nlinks[i]) {
								emembership[e] = nlinks[i] - eBegin;
								break;
							}
						} else {
							counter = 1;
						}
					}
				}
			}

			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				std::vector<eslocal>  tboundaryEDistribution;
				std::vector<eslocal>  tboundaryEData;
				std::vector<Element*> tboundaryEPointers;
				if (t == 0) {
					tboundaryEDistribution.push_back(0);
				}

				std::vector<size_t> tsize(targets.size());
				for (size_t i = 0; i < targets.size(); i++) {
					tsize[i] = sBoundary[t][i].size();
				}
				eslocal target;
				auto enodes = _mesh->boundaryRegions[r]->elements->cbegin() + distribution[t];
				const auto &IDs = _mesh->nodes->IDs->datatarray();
				const auto &epointer = _mesh->boundaryRegions[r]->epointers->datatarray();
				for (size_t e = distribution[t]; e < distribution[t + 1]; ++e, ++enodes) {
					if (partition[emembership[e]] == environment->MPIrank) {
						tboundaryEPointers.push_back(epointer[e]);
						for (auto n = enodes->begin(); n != enodes->end(); ++n) {
							tboundaryEData.push_back(IDs[*n]);
						}
						tboundaryEDistribution.push_back(tboundaryEData.size());
					} else {
						target = t2i(partition[emembership[e]]);
						sBoundary[t][target].push_back(static_cast<int>(epointer[e]->code));
						for (auto n = enodes->begin(); n != enodes->end(); ++n) {
							sBoundary[t][target].push_back(IDs[*n]);
						}
					}
				}
				for (size_t i = 0; i < targets.size(); i++) {
					sBoundary[t][i][r] = sBoundary[t][i].size() - tsize[i];
				}

				boundaryEDistribution[r][t].swap(tboundaryEDistribution);
				boundaryEData[r][t].swap(tboundaryEData);
				boundaryEPointers[r][t].swap(tboundaryEPointers);
			}
		}
	}

	e21.end();
	eval.addEvent(e21);

	TimeEvent e3("send data");
	e3.start();

	// Step 3: Send data to target processes

	#pragma omp parallel for
	for (size_t target = 0; target < targets.size(); ++target) {
		sElements[0][target].front() = sElements[0][target].size();
		sNodes[0][target].front() = sNodes[0][target].size();
		for (size_t t = 1; t < threads; t++) {
			sElements[t][target].front() = sElements[t][target].size();
			sElements[0][target].insert(sElements[0][target].end(), sElements[t][target].begin(), sElements[t][target].end());
			sNodes[t][target].front() = sNodes[t][target].size();
			sNodes[0][target].insert(sNodes[0][target].end(), sNodes[t][target].begin(), sNodes[t][target].end());
			sBoundary[0][target].insert(sBoundary[0][target].end(), sBoundary[t][target].begin(), sBoundary[t][target].end());
		}
	}

	if (!Communication::sendVariousTargets(sElements[0], rElements, targets)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange elements data.";
	}

	if (!Communication::sendVariousTargets(sNodes[0], rNodes, targets)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange nodes data.";
	}

	if (!Communication::sendVariousTargets(sBoundary[0], rBoundary, targets)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange boundary data.";
	}

	e3.end();
	eval.addEvent(e3);

	TimeEvent e4("deserialize elements");
	e4.start();

	// Step 4: Deserialize element data
	for (size_t i = 0; i < rElements.size(); i++) {
		std::vector<size_t> rdistribution({ 0 });
		size_t p = 0;
		while (p < rElements[i].size()) {
			rdistribution.push_back(p += rElements[i][p]);
		}

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {

			std::vector<eslocal>  telemsIDs;
			std::vector<int>      telemsBody;
			std::vector<int>      telemsMaterial;
			std::vector<Element*> telemsEpointer;
			std::vector<eslocal>  telemsNodesDistribution;
			std::vector<eslocal>  telemsNodesData;
			std::vector<eslocal>  telemsRegions;

			telemsIDs.reserve(rdistribution[t + 1] - rdistribution[t]);
			telemsBody.reserve(rdistribution[t + 1] - rdistribution[t]);
			telemsMaterial.reserve(rdistribution[t + 1] - rdistribution[t]);
			telemsEpointer.reserve(rdistribution[t + 1] - rdistribution[t]);
			telemsNodesDistribution.reserve(rdistribution[t + 1] - rdistribution[t] + 1);
			telemsRegions.reserve(rdistribution[t + 1] - rdistribution[t]);

			eslocal distOffset = 0;

			if (elemsNodesDistribution[t].size()) {
				distOffset = elemsNodesDistribution[t].back();
			}
			if (t == 0 && elemsNodesDistribution[t].size() == 0) {
				telemsNodesDistribution.push_back(0);
			}

			for (size_t e = rdistribution[t] + 1; e < rdistribution[t + 1]; ) {
				telemsIDs.push_back(rElements[i][e++]);
				telemsBody.push_back(rElements[i][e++]);
				telemsMaterial.push_back(rElements[i][e++]);
				telemsEpointer.push_back(_mesh->_eclasses[t] + rElements[i][e++]);
				telemsNodesData.insert(telemsNodesData.end(), rElements[i].begin() + e + 1, rElements[i].begin() + e + 1 + rElements[i][e]);
				telemsNodesDistribution.push_back(telemsNodesData.size() + distOffset);
				e += rElements[i][e++]; // nodes + nodes size
				telemsRegions.insert(telemsRegions.end(), rElements[i].begin() + e, rElements[i].begin() + e + eregionsBitMaskSize);
				e += eregionsBitMaskSize;
			}

			elemsIDs[t].insert(elemsIDs[t].end(), telemsIDs.begin(), telemsIDs.end());
			elemsBody[t].insert(elemsBody[t].end(), telemsBody.begin(), telemsBody.end());
			elemsMaterial[t].insert(elemsMaterial[t].end(), telemsMaterial.begin(), telemsMaterial.end());
			elemsEpointer[t].insert(elemsEpointer[t].end(), telemsEpointer.begin(), telemsEpointer.end());
			elemsNodesDistribution[t].insert(elemsNodesDistribution[t].end(), telemsNodesDistribution.begin(), telemsNodesDistribution.end());
			elemsNodesData[t].insert(elemsNodesData[t].end(), telemsNodesData.begin(), telemsNodesData.end());
			elemsRegions[t].insert(elemsRegions[t].end(), telemsRegions.begin(), telemsRegions.end());
		}
	}

	e4.end();
	eval.addEvent(e4);

	TimeEvent e41("deserialize nodes");
	e41.start();

	// Step 4: Deserialize node data
	std::vector<eslocal> nodeset;
	for (size_t t = 0; t < threads; t++) {
		nodeset.insert(nodeset.end(), nodesIDs[t].begin(), nodesIDs[t].end());
	}
	std::sort(nodeset.begin(), nodeset.end());
	for (size_t i = 0; i < rNodes.size(); i++) {
		std::vector<size_t> rdistribution({ 0 });
		size_t p = 0;
		while (p < rNodes[i].size()) {
			rdistribution.push_back(p += rNodes[i][p]);
		}
		std::vector<std::vector<eslocal> > tnodeset(threads);
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			Point point;

			std::vector<eslocal>  tnodesIDs;
			std::vector<Point>    tnodesCoordinates;
			std::vector<eslocal>  tnodesElemsDistribution;
			std::vector<eslocal>  tnodesElemsData;
			std::vector<eslocal>  tnodesRegions;
			std::vector<eslocal>  tnodeSet;

			tnodesIDs.reserve(rdistribution[t + 1] - rdistribution[t]);
			tnodesCoordinates.reserve(rdistribution[t + 1] - rdistribution[t]);
			tnodesElemsDistribution.reserve(rdistribution[t + 1] - rdistribution[t]);
			tnodesRegions.reserve(rdistribution[t + 1] - rdistribution[t]);
			tnodeSet.reserve(rdistribution[t + 1] - rdistribution[t]);

			eslocal distOffset = 0;

			if (nodesElemsDistribution[t].size()) {
				distOffset = nodesElemsDistribution[t].back();
			}
			if (t == 0 && nodesElemsDistribution[t].size() == 0) {
				tnodesElemsDistribution.push_back(0);
			}

			auto nodesetit = nodeset.begin();
			if (rdistribution[t] + 1 < rdistribution[t + 1]) {
				nodesetit = std::lower_bound(nodeset.begin(), nodeset.end(), rNodes[i][rdistribution[t] + 1]);
			}
			for (size_t n = rdistribution[t] + 1; n < rdistribution[t + 1]; ) {
				while (nodesetit != nodeset.end() && *nodesetit < rNodes[i][n]) ++nodesetit;
				if (nodesetit == nodeset.end() || *nodesetit != rNodes[i][n]) {
					tnodesIDs.push_back(rNodes[i][n]);
					tnodeSet.push_back(rNodes[i][n]);
					n += 1; //ID
					memcpy(&point, rNodes[i].data() + n, sizeof(Point));
					tnodesCoordinates.push_back(point);
					n += sizeof(Point) / sizeof(eslocal); // points
					tnodesElemsData.insert(tnodesElemsData.end(), rNodes[i].begin() + n + 1, rNodes[i].begin() + n + 1 + rNodes[i][n]);
					tnodesElemsDistribution.push_back(tnodesElemsData.size() + distOffset);
					n += rNodes[i][n] + 1; // linksize + links
					tnodesRegions.insert(tnodesRegions.end(), rNodes[i].begin() + n, rNodes[i].begin() + n + bregionsBitMaskSize);
					n += bregionsBitMaskSize; // region mask
				} else {
					n += 1 + sizeof(Point) / sizeof(eslocal); // id, Point
					n += 1 + rNodes[i][n]; // linksize, links
					n += bregionsBitMaskSize; // region mask
				}
			}

			nodesIDs[t].insert(nodesIDs[t].end(), tnodesIDs.begin(), tnodesIDs.end());
			nodesCoordinates[t].insert(nodesCoordinates[t].end(), tnodesCoordinates.begin(), tnodesCoordinates.end());
			nodesElemsDistribution[t].insert(nodesElemsDistribution[t].end(), tnodesElemsDistribution.begin(), tnodesElemsDistribution.end());
			nodesElemsData[t].insert(nodesElemsData[t].end(), tnodesElemsData.begin(), tnodesElemsData.end());
			nodesRegions[t].insert(nodesRegions[t].end(), tnodesRegions.begin(), tnodesRegions.end());
			tnodeset[t].swap(tnodeSet);
		}

		size_t nsize = nodeset.size();
		for (size_t t = 0; t < threads; t++) {
			nodeset.insert(nodeset.end(), tnodeset[t].begin(), tnodeset[t].end());
		}
		std::inplace_merge(nodeset.begin(), nodeset.begin() + nsize, nodeset.end());
	}

	Esutils::threadDistributionToFullDistribution(elemsNodesDistribution);
	Esutils::threadDistributionToFullDistribution(nodesElemsDistribution);

	e41.end();
	eval.addEvent(e41);

	TimeEvent e42("deserialize boundary");
	e42.start();

	// Step 4: Deserialize boundary data
	for (size_t n = 0; n < rBoundary.size(); ++n) {
		std::vector<std::vector<eslocal> > toffset(_mesh->boundaryRegions.size()), tsize(_mesh->boundaryRegions.size());
		eslocal offset = 0, p = 0;
		for (size_t t = 0; t < threads; t++) {
			for (size_t r = 0; r < _mesh->boundaryRegions.size(); r++) {
				toffset[r].push_back(offset + _mesh->boundaryRegions.size());
				tsize[r].push_back(rBoundary[n][p + r]);
				offset += rBoundary[n][p + r];
			}
			offset += _mesh->boundaryRegions.size();
			p = offset;
		}
		for (size_t r = 0; r < _mesh->boundaryRegions.size(); r++) {
			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {

				std::vector<eslocal>  tboundaryEDistribution;
				std::vector<eslocal>  tboundaryEData;
				std::vector<Element*> tboundaryEPointers;

				eslocal distOffset = 0;

				if (boundaryEDistribution[r][t].size()) {
					distOffset = boundaryEDistribution[r][t].back();
				}
				if (t == 0 && boundaryEDistribution[r][t].size() == 0) {
					tboundaryEDistribution.push_back(0);
				}

				for (eslocal i = toffset[r][t]; i < toffset[r][t] + tsize[r][t];) {
					tboundaryEPointers.push_back(&_mesh->_eclasses[t][rBoundary[n][i++]]);
					tboundaryEData.insert(tboundaryEData.end(), rBoundary[n].begin() + i, rBoundary[n].begin() + i + tboundaryEPointers.back()->nodes);
					tboundaryEDistribution.push_back(tboundaryEData.size() + distOffset);
					i += tboundaryEPointers.back()->nodes;
				}

				boundaryEPointers[r][t].insert(boundaryEPointers[r][t].end(), tboundaryEPointers.begin(), tboundaryEPointers.end());
				boundaryEDistribution[r][t].insert(boundaryEDistribution[r][t].end(), tboundaryEDistribution.begin(), tboundaryEDistribution.end());
				boundaryEData[r][t].insert(boundaryEData[r][t].end(), tboundaryEData.begin(), tboundaryEData.end());
			}
		}
	}

	// elements are redistributed later while decomposition -> distribution is not changed now
	std::vector<size_t> elemDistribution(threads + 1);
	for (size_t t = 1; t <= threads; t++) {
		elemDistribution[t] = elemDistribution[t - 1] + elemsIDs[t - 1].size();
	}

	elements->IDs = new serializededata<eslocal, eslocal>(1, elemsIDs);
	elements->body = new serializededata<eslocal, int>(1, elemsBody);
	elements->material = new serializededata<eslocal, int>(1, elemsMaterial);
	elements->epointers = new serializededata<eslocal, Element*>(1, elemsEpointer);
	elements->nodes = new serializededata<eslocal, eslocal>(elemsNodesDistribution, elemsNodesData); // global IDs

	elements->regionMaskSize = eregionsBitMaskSize;
	elements->regions = new serializededata<eslocal, int>(eregionsBitMaskSize, elemsRegions);

	elements->size = elements->IDs->structures();
	elements->distribution = elements->IDs->datatarray().distribution();

	// Step 5: Balance node data to threads
	std::vector<size_t> nodeDistribution(threads);
	for (size_t t = 1; t < threads; t++) {
		nodeDistribution[t] = nodeDistribution[t - 1] + nodesIDs[t - 1].size();
	}

	serializededata<eslocal, eslocal>::balance(1, nodesIDs);
	serializededata<eslocal, Point>::balance(1, nodesCoordinates);
	serializededata<eslocal, eslocal>::balance(nodesElemsDistribution, nodesElemsData);

	nodes->IDs = new serializededata<eslocal, eslocal>(1, nodesIDs);
	nodes->coordinates = new serializededata<eslocal, Point>(1, nodesCoordinates);
	nodes->elements = new serializededata<eslocal, eslocal>(nodesElemsDistribution, nodesElemsData);
	nodes->size = nodes->IDs->datatarray().size();
	nodes->distribution = nodes->IDs->datatarray().distribution();

	e42.end();
	eval.addEvent(e42);

	TimeEvent e5("balance nodes");
	e5.start();

	for (size_t r = 0; r < _mesh->boundaryRegions.size(); r++) {
		if (_mesh->boundaryRegions[r]->dimension) {
			delete _mesh->boundaryRegions[r]->elements;
			delete _mesh->boundaryRegions[r]->epointers;

			Esutils::threadDistributionToFullDistribution(boundaryEDistribution[r]);
			_mesh->boundaryRegions[r]->elements = new serializededata<eslocal, eslocal>(boundaryEDistribution[r], boundaryEData[r]);
			_mesh->boundaryRegions[r]->epointers = new serializededata<eslocal, Element*>(1, boundaryEPointers[r]);
		}
	}

	for (size_t r = 0; r < _mesh->elementsRegions.size(); r++) {
		int maskOffset = r / (8 * sizeof(eslocal));
		int bit = 1 << (r % (8 * sizeof(eslocal)));
		delete _mesh->elementsRegions[r]->elements;
		std::vector<std::vector<eslocal> > regionelems(threads);

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (size_t i = 0; i < elemsRegions[t].size(); i += eregionsBitMaskSize) {
				if (elemsRegions[t][i + maskOffset] & bit) {
					regionelems[t].push_back(elemDistribution[t] + i / eregionsBitMaskSize);
				}
			}
		}

		serializededata<eslocal, eslocal>::balance(1, regionelems);
		_mesh->elementsRegions[r]->elements = new serializededata<eslocal, eslocal>(1, regionelems);
	}

	for (size_t r = 0; r < _mesh->boundaryRegions.size(); r++) {
		if (_mesh->boundaryRegions[r]->nodes) {
			int maskOffset = r / (8 * sizeof(eslocal));
			int bit = 1 << (r % (8 * sizeof(eslocal)));
			delete _mesh->boundaryRegions[r]->nodes;
			std::vector<std::vector<eslocal> > regionnodes(threads);

			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				for (size_t i = 0; i < nodesRegions[t].size(); i += bregionsBitMaskSize) {
					if (nodesRegions[t][i + maskOffset] & bit) {
						regionnodes[t].push_back(nodeDistribution[t] + i / bregionsBitMaskSize);
					}
				}
			}

			serializededata<eslocal, eslocal>::balance(1, regionnodes);
			_mesh->boundaryRegions[r]->nodes = new serializededata<eslocal, eslocal>(1, regionnodes);
		}
	}

	e5.end();
	eval.addEvent(e5);

	TimeEvent e__("NEW REINDEXING");
	e__.start();

	std::vector<eslocal> eIDsOLD = _mesh->elements->gatherElementsProcDistribution();
	std::vector<eslocal> eIDsNEW = elements->gatherElementsProcDistribution();

	for (size_t t = 1; t < threads; ++t) {
		elemsIDs[0].insert(elemsIDs[0].end(), elemsIDs[t].begin(), elemsIDs[t].end());
	}

	std::vector<eslocal> epermutation(elemsIDs[0].size());
	std::iota(epermutation.begin(), epermutation.end(), 0);
	std::sort(epermutation.begin(), epermutation.end(), [&] (eslocal i, eslocal j) { return elemsIDs[0][i] < elemsIDs[0][j]; });

	Esutils::sortWithInplaceMerge(elemsIDs[0], elemDistribution);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; ++t) {
		Esutils::sortAndRemoveDuplicity(nodesElemsData[t]);
	}
	Esutils::inplaceMerge(nodesElemsData);
	Esutils::sortAndRemoveDuplicity(nodesElemsData);

	std::vector<std::vector<eslocal> > requestedIDs, receivedTargets, IDrequests;
	std::vector<std::vector<int> > IDtargets(threads);
	std::vector<int> sources;

	std::vector<size_t> rdistribution = tarray<int>::distribute(threads, environment->MPIsize);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; ++t) {
		std::vector<int> ttargets;

		auto eIDbegin = std::lower_bound(nodesElemsData[0].begin(), nodesElemsData[0].end(), eIDsOLD[rdistribution[t]]);
		auto eIDend   = std::lower_bound(eIDbegin, nodesElemsData[0].end(), eIDsOLD[rdistribution[t + 1]]);

		for (size_t r = rdistribution[t]; r < rdistribution[t + 1]; r++) {
			eIDbegin = std::lower_bound(eIDbegin, eIDend, eIDsOLD[r]);
			auto end = std::lower_bound(eIDbegin, eIDend, eIDsOLD[r + 1]);
			if (eIDbegin != end) {
				ttargets.push_back(r);
			}
			eIDbegin = end;
		}

		IDtargets[t].swap(ttargets);
	}

	for (size_t t = 1; t < threads; ++t) {
		IDtargets[0].insert(IDtargets[0].end(), IDtargets[t].begin(), IDtargets[t].end());
	}

	requestedIDs.resize(IDtargets[0].size());

	#pragma omp parallel for
	for (size_t t = 0; t < IDtargets[0].size(); t++) {
		auto mybegin = std::lower_bound(elemsIDs[0].begin(), elemsIDs[0].end(), eIDsOLD[IDtargets[0][t]]);
		auto myend   = std::lower_bound(elemsIDs[0].begin(), elemsIDs[0].end(), eIDsOLD[IDtargets[0][t] + 1]);
		auto nbegin = std::lower_bound(nodesElemsData[0].begin(), nodesElemsData[0].end(), eIDsOLD[IDtargets[0][t]]);
		auto nend   = std::lower_bound(nodesElemsData[0].begin(), nodesElemsData[0].end(), eIDsOLD[IDtargets[0][t] + 1]);
		requestedIDs[t].resize(nend - nbegin);
		requestedIDs[t].resize(std::set_difference(nbegin, nend, mybegin, myend, requestedIDs[t].begin()) - requestedIDs[t].begin());
	}

	if (!Communication::sendVariousTargets(requestedIDs, IDrequests, IDtargets[0], sources)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange ID requests.";
	}

	for (size_t r = 0; r < IDrequests.size(); r++) {
		rdistribution = tarray<size_t>::distribute(threads, IDrequests[r].size());

		#pragma omp parallel for
		for (size_t t = 0; t < threads; ++t) {
			if (rdistribution[t] != rdistribution[t + 1]) {
				for (size_t e = rdistribution[t]; e < rdistribution[t + 1]; ++e) {
					IDrequests[r][e] = partition[IDrequests[r][e] - eIDsOLD[environment->MPIrank]];
				}
			}
		}
	}

	if (!Communication::sendVariousTargets(IDrequests, receivedTargets, sources)) {
		ESINFO(ERROR) << "ESPRESO internal error: return ID targets.";
	}

	IDtargets.clear();
	IDtargets.resize(threads);
	for (size_t r = 0; r < receivedTargets.size(); r++) {
		rdistribution = tarray<size_t>::distribute(threads, receivedTargets[r].size());

		#pragma omp parallel for
		for (size_t t = 0; t < threads; ++t) {
			std::vector<int> ttargets;

			std::vector<int>::iterator it;
			for (size_t e = rdistribution[t]; e < rdistribution[t + 1]; ++e) {
				if (receivedTargets[r][e] != environment->MPIrank) {
					it = std::lower_bound(ttargets.begin(), ttargets.end(), receivedTargets[r][e]);
					if (it == ttargets.end() || *it != receivedTargets[r][e]) {
						ttargets.insert(it, receivedTargets[r][e]);
					}
				}
			}

			IDtargets[t].insert(IDtargets[t].end(), ttargets.begin(), ttargets.end());
		}
	}
	for (size_t t = 0; t < threads; ++t) {
		Esutils::sortAndRemoveDuplicity(IDtargets[t]);
	}
	for (size_t t = 1; t < threads; ++t) {
		IDtargets[0].insert(IDtargets[0].end(), IDtargets[t].begin(), IDtargets[t].end());
	}
	Esutils::sortAndRemoveDuplicity(IDtargets[0]);

	std::vector<std::vector<std::vector<eslocal> > > newIDrequests(threads, std::vector<std::vector<eslocal> >(IDtargets[0].size()));
	std::vector<std::vector<eslocal> > newIDs;

	for (size_t r = 0; r < receivedTargets.size(); r++) {
		rdistribution = tarray<size_t>::distribute(threads, receivedTargets[r].size());

		#pragma omp parallel for
		for (size_t t = 0; t < threads; ++t) {
			size_t tindex;
			std::vector<std::vector<eslocal> > tnewIDrequests(IDtargets[0].size());

			for (size_t e = rdistribution[t]; e < rdistribution[t + 1]; ++e) {
				tindex = std::lower_bound(IDtargets[0].begin(), IDtargets[0].end(), receivedTargets[r][e]) - IDtargets[0].begin();
				tnewIDrequests[tindex].push_back(requestedIDs[r][e]);
			}

			for (size_t n = 0; n < IDtargets[0].size(); n++) {
				newIDrequests[t][n].insert(newIDrequests[t][n].end(), tnewIDrequests[n].begin(), tnewIDrequests[n].end());
			}
		}

		for (size_t t = 1; t < threads; ++t) {
			for (size_t n = 0; n < IDtargets[0].size(); n++) {
				newIDrequests[0][n].insert(newIDrequests[0][n].end(), newIDrequests[t][n].begin(), newIDrequests[t][n].end());
				newIDrequests[t][n].clear();
			}
		}
	}


	if (!Communication::sendVariousTargets(newIDrequests[0], IDrequests, IDtargets[0], sources)) {
		ESINFO(ERROR) << "ESPRESO internal error: request new ID.";
	}

	for (size_t r = 0; r < IDrequests.size(); r++) {
		rdistribution = tarray<size_t>::distribute(threads, IDrequests[r].size());

		#pragma omp parallel for
		for (size_t t = 0; t < threads; ++t) {
			if (rdistribution[t] != rdistribution[t + 1]) {
				auto eit = std::lower_bound(elemsIDs[0].begin(), elemsIDs[0].end(), IDrequests[r][rdistribution[t]]);
				for (size_t e = rdistribution[t]; e < rdistribution[t + 1]; ++e) {
					while (eit != elemsIDs[0].end() && *eit < IDrequests[r][e]) ++eit;
					IDrequests[r][e] = epermutation[eit - elemsIDs[0].begin()] + eIDsNEW[environment->MPIrank];
				}
			}
		}
	}

	if (!Communication::sendVariousTargets(IDrequests, newIDs, sources)) {
		ESINFO(ERROR) << "ESPRESO internal error: return new ID.";
	}

	std::vector<size_t> offsets = { 0 };
	for (size_t r = 0; r < newIDs.size(); r++) {
		offsets.push_back(offsets.back() + newIDs[r].size());
	}
	std::vector<std::pair<eslocal, eslocal> > IDMap(offsets.back());

	for (size_t r = 0; r < newIDs.size(); r++) {
		rdistribution = tarray<size_t>::distribute(threads, newIDrequests[0][r].size());

		#pragma omp parallel for
		for (size_t t = 0; t < threads; ++t) {
			for (size_t e = rdistribution[t]; e < rdistribution[t + 1]; ++e) {
				IDMap[offsets[r] + e].first = newIDrequests[0][r][e];
				IDMap[offsets[r] + e].second = newIDs[r][e];
			}
		}
	}

	Esutils::sortWithInplaceMerge(IDMap, offsets);

	e__.end();
	eval.addEvent(e__);

	TimeEvent e15("reindex node->elements");
	e15.start();

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (auto elem = nodes->elements->begin(t); elem != nodes->elements->end(t); ++elem) {
			for (auto e = elem->begin(); e != elem->end(); ++e) {
				auto it = std::lower_bound(IDMap.begin(), IDMap.end(), std::make_pair(*e, 0));
				if (it == IDMap.end() || it->first != *e) {
					*e = eIDsNEW[environment->MPIrank] + epermutation[std::lower_bound(elemsIDs[0].begin(), elemsIDs[0].end(), *e) - elemsIDs[0].begin()];
				} else {
					*e = it->second;
				}
			}
			std::sort(elem->begin(), elem->end());
		}
	}

	e15.end();
	eval.addEvent(e15);

	TimeEvent e16("reindex elements->nodes");
	e16.start();

	std::vector<eslocal> permutation(nodes->size);
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) { return nodes->IDs->datatarray().data()[i] < nodes->IDs->datatarray().data()[j]; });

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (auto n = elements->nodes->begin(t)->begin(); n != elements->nodes->end(t)->begin(); ++n) {
			*n = *std::lower_bound(permutation.begin(), permutation.end(), *n, [&] (eslocal i, eslocal val) {
				return nodes->IDs->datatarray().data()[i] < val;
			});
		}
	}

	for (size_t r = 0; r < _mesh->boundaryRegions.size(); r++) {
		if (_mesh->boundaryRegions[r]->elements) {
			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				for (auto n = _mesh->boundaryRegions[r]->elements->begin(t)->begin(); n != _mesh->boundaryRegions[r]->elements->end(t)->begin(); ++n) {
					*n = *std::lower_bound(permutation.begin(), permutation.end(), *n, [&] (eslocal i, eslocal val) {
						return nodes->IDs->datatarray().data()[i] < val;
					});
				}
			}
		}
	}

	e16.end();
	eval.addEvent(e16);

	TimeEvent e17("reindex ranks");
	e17.start();

	std::vector<std::vector<eslocal> > rankBoundaries(threads);
	std::vector<std::vector<int> > rankData(threads);

	rankBoundaries.front().push_back(0);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		int rank, rsize;
		for (auto elem = nodes->elements->begin(t); elem != nodes->elements->end(t); ++elem) {
			rsize = 1;
			rankData[t].push_back(std::lower_bound(eIDsNEW.begin(), eIDsNEW.end(), *elem->begin() + 1) - eIDsNEW.begin() - 1);
			for (auto e = elem->begin() + 1; e != elem->end(); ++e) {
				rank = std::lower_bound(eIDsNEW.begin(), eIDsNEW.end(), *e + 1) - eIDsNEW.begin() - 1;
				if (rank != rankData[t].back()) {
					rankData[t].push_back(rank);
					++rsize;
				}
			}
			rankBoundaries[t].push_back(rsize);
			if (rankBoundaries[t].size() > 1) {
				rankBoundaries[t].back() += *(rankBoundaries[t].end() - 2);
			}
		}
	}

	Esutils::threadDistributionToFullDistribution(rankBoundaries);

	nodes->ranks = new serializededata<eslocal, int>(rankBoundaries, rankData);

	e17.end();
	eval.addEvent(e17);

	TimeEvent e18("swap");
	e18.start();

	std::iota(elements->IDs->datatarray().begin(), elements->IDs->datatarray().end(), eIDsNEW[environment->MPIrank]);
	std::swap(_mesh->elements, elements);
	std::swap(_mesh->nodes, nodes);
	delete _mesh->halo;
	_mesh->halo = new ElementStore(_mesh->_eclasses);
	_mesh->neighbours.clear();
	for (size_t t = 0; t < IDtargets[0].size(); t++) {
		if (IDtargets[0][t] != environment->MPIrank) {
			_mesh->neighbours.push_back(IDtargets[0][t]);
		}
	}
	_mesh->neighboursWithMe = _mesh->neighbours;
	_mesh->neighboursWithMe.push_back(environment->MPIrank);
	std::sort(_mesh->neighboursWithMe.begin(), _mesh->neighboursWithMe.end());

	delete elements;
	delete nodes;

	e18.end();
	eval.addEvent(e18);

//	eval.totalTime.endWithBarrier();
//	eval.printStatsMPI();
//	exit(0);

	this->computeElementsNeighbors();

	finish("exchanging elements");
}

void MeshPreprocessing::permuteElements(const std::vector<eslocal> &permutation, const std::vector<size_t> &distribution)
{
	start("permutation of elements");

	if (_mesh->nodes->elements == NULL) {
		this->linkNodesAndElements();
	}

	std::vector<eslocal> backpermutation(permutation.size());
	std::iota(backpermutation.begin(), backpermutation.end(), 0);
	std::sort(backpermutation.begin(), backpermutation.end(), [&] (eslocal i, eslocal j) { return permutation[i] < permutation[j]; });

	size_t threads = environment->OMP_NUM_THREADS;

	auto n2i = [ & ] (size_t neighbor) {
		return std::lower_bound(_mesh->neighbours.begin(), _mesh->neighbours.end(), neighbor) - _mesh->neighbours.begin();
	};

	std::vector<eslocal> IDBoundaries = _mesh->elements->gatherElementsProcDistribution();
	std::vector<std::vector<std::pair<eslocal, eslocal> > > rHalo(_mesh->neighbours.size());

	if (_mesh->elements->neighbors != NULL || _mesh->nodes->elements != NULL) {
		// thread x neighbor x elements(oldID, newID)
		std::vector<std::vector<std::vector<std::pair<eslocal, eslocal> > > > sHalo(threads, std::vector<std::vector<std::pair<eslocal, eslocal> > >(_mesh->neighbours.size()));

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			auto ranks = _mesh->nodes->ranks->cbegin(t);
			auto elements = _mesh->nodes->elements->cbegin(t);
			eslocal begine = IDBoundaries[environment->MPIrank];
			eslocal ende   = IDBoundaries[environment->MPIrank + 1];

			for (auto n = _mesh->nodes->distribution[t]; n < _mesh->nodes->distribution[t + 1]; ++n, ++ranks, ++elements) {
				for (auto rank = ranks->begin(); rank != ranks->end(); ++rank) {
					if (*rank != environment->MPIrank) {
						for (auto e = elements->begin(); e != elements->end(); ++e) {
							if (begine <= *e && *e < ende) {
								sHalo[t][n2i(*rank)].push_back(std::make_pair(*e, backpermutation[*e - begine] + begine));
							}
						}
					}
				}
			}

			for (size_t n = 0; n < sHalo[t].size(); ++n) {
				Esutils::sortAndRemoveDuplicity(sHalo[t][n]);
			}
		}

		Esutils::mergeThreadedUniqueData(sHalo);

		if (!Communication::exchangeUnknownSize(sHalo[0], rHalo, _mesh->neighbours)) {
			ESINFO(ERROR) << "ESPRESO internal error: exchange halo element new IDs while element permutation.";
		}
	}

	auto globalremap = [&] (serializededata<eslocal, eslocal>* data, bool sort) {
		if (data == NULL) {
			return;
		}
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			int source;
			for (auto e = data->begin(t); e != data->end(t); ++e) {
				for (auto n = e->begin(); n != e->end(); ++n) {
					if (*n >= 0) {
						source = std::lower_bound(IDBoundaries.begin(), IDBoundaries.end(), *n + 1) - IDBoundaries.begin() - 1;
						if (source == environment->MPIrank) {
							*n = IDBoundaries[environment->MPIrank] + backpermutation[*n - IDBoundaries[environment->MPIrank]];
						} else {
							*n = std::lower_bound(rHalo[n2i(source)].begin(), rHalo[n2i(source)].end(), std::make_pair(*n, 0))->second;
						}
					}
				}
				if (sort) {
					std::sort(e->begin(), e->end());
				}
			}
		}
	};

	eslocal firstID = _mesh->elements->IDs->datatarray().front();
	_mesh->elements->permute(permutation, distribution);
	std::iota(_mesh->elements->IDs->datatarray().begin(), _mesh->elements->IDs->datatarray().end(), firstID);

	globalremap(_mesh->elements->neighbors, false);
	globalremap(_mesh->nodes->elements, true);

	for (size_t r = 0; r < _mesh->elementsRegions.size(); ++r) {
		for (auto n = _mesh->elementsRegions[r]->elements->datatarray().begin(); n != _mesh->elementsRegions[r]->elements->datatarray().end(); ++n) {
			*n = backpermutation[*n];
		}
		std::sort(_mesh->elementsRegions[r]->elements->datatarray().begin(), _mesh->elementsRegions[r]->elements->datatarray().end());
	}

	finish("permutation of elements");
}
