
#include "meshpreprocessing.h"

#include "mesh/mesh.h"

#include "mesh/store/store.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementsregionstore.h"
#include "mesh/store/boundaryregionstore.h"

#include "mesh/elements/element.h"

#include "globals/run.h"

#include "basis/containers/point.h"
#include "basis/containers/serializededata.h"
#include "basis/utilities/communication.h"
#include "basis/utilities/utils.h"
#include "basis/utilities/parser.h"
#include "basis/logging/logging.h"
#include "basis/logging/timeeval.h"

#include "config/ecf/root.h"
#include "config/ecf/decomposition.h"

#include <algorithm>
#include <numeric>
#include <cstring>

#include "wrappers/metis/metiswrapper.h"
#include "wrappers/metis/parmetiswrapper.h"

using namespace espreso;


void MeshPreprocessing::reclusterize()
{
	if (environment->MPIsize == 1) {
		skip("re-distribution of the mesh to processes");
		return;
	}

	if (_mesh->elements->neighbors == NULL) {
		this->computeElementsNeighbors();
	}

	if (_mesh->halo->IDs == NULL) {
		exchangeHalo();
	}

	// Disable due to horible scalability
//	if (_mesh->elements->centers == NULL) {
//		computeElementsCenters();
//	}

	start("compute global dual graph");

	bool separateRegions = run::ecf->decomposition.separate_regions;
	bool separateMaterials = run::ecf->decomposition.separate_materials;
	bool separateEtypes = run::ecf->decomposition.separate_etypes;

	if (separateRegions && _mesh->elements->regions == NULL) {
		fillRegionMask();
	}

	size_t threads = environment->OMP_NUM_THREADS;

	esint eoffset = _mesh->elements->gatherElementsProcDistribution()[environment->MPIrank];

	std::vector<esint> dDistribution(_mesh->elements->size + 1);
	std::vector<std::vector<esint> > dData(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<esint> tdata;
		int mat1 = 0, mat2 = 0, reg = 0, etype1 = 0, etype2 = 0;
		int rsize = _mesh->elements->regionMaskSize;
		esint hindex;

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
							reg = memcmp(_mesh->elements->regions->datatarray().data() + e * rsize, _mesh->halo->regions->datatarray().data() + hindex * rsize, sizeof(esint) * rsize);
						} else {
							reg = memcmp(_mesh->elements->regions->datatarray().data() + e * rsize, _mesh->elements->regions->datatarray().data() + (*n - eoffset) * rsize, sizeof(esint) * rsize);
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

	std::vector<esint> partition(_mesh->elements->size, environment->MPIrank);

	finish("compute global dual graph");

	MPISubset subset(run::ecf->decomposition.metis_options, MPITools::procs());

	start("ParMETIS::KWay");
	esint edgecut = ParMETIS::call(
			ParMETIS::METHOD::ParMETIS_V3_PartKway, subset,
			dDistribution, dData.front(), partition
	);
	finish("ParMETIS::KWay");

	if (run::ecf->decomposition.metis_options.refinement) {
		start("ParMETIS::AdaptiveRepart");
		esint prev = 2 * edgecut;
		while (1.01 * edgecut < prev) {
			prev = edgecut;
			edgecut = ParMETIS::call(
					ParMETIS::METHOD::ParMETIS_V3_RefineKway, subset,
					dDistribution, dData.front(), partition
			);
		}
		finish("ParMETIS::AdaptiveRepart");
	}

	this->exchangeElements(partition);
}

void MeshPreprocessing::partitiate(esint parts)
{
	std::vector<esint> dualDist, dualData;
	this->computeDecomposedDual(dualDist, dualData);

	start("check continuity");

	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<int> partID(_mesh->elements->size, -1);

	int nextID = 0;
	for (esint e = 0; e < _mesh->elements->size; ++e) {
		std::vector<esint> stack;
		if (partID[e] == -1) {
			stack.push_back(e);
			partID[e] = nextID;
			while (stack.size()) {
				esint current = stack.back();
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
	std::vector<esint> partition(_mesh->elements->size);

	finish("check continuity");

	if (nextID == 1) {
		start("process non-continuous dual graph");
		finish("process non-continuous dual graph");

		start("METIS::KWay");
		METIS::call(
				run::ecf->decomposition.metis_options,
				_mesh->elements->size,
				dualDist.data(), dualData.data(),
				0, NULL, NULL,
				parts, partition.data());
		finish("METIS::KWay");
		clusters.resize(parts, 0);
		_mesh->elements->nclusters = 1;

		start("reindex METIS output");
		finish("reindex METIS output");

	} else { // non-continuous dual graph
		start("process non-continuous dual graph");

		// thread x part x elements
		std::vector<std::vector<std::vector<esint> > > tdecomposition(threads, std::vector<std::vector<esint> >(nextID));
		std::vector<std::vector<esint> > tdualsize(threads, std::vector<esint>(nextID));

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (size_t e = _mesh->elements->distribution[t]; e < _mesh->elements->distribution[t + 1]; ++e) {
				tdecomposition[t][partID[e]].push_back(e);
				tdualsize[t][partID[e]] += dualDist[e + 1] - dualDist[e];
			}
		}
		std::vector<std::vector<esint> > foffsets(nextID), noffsets(nextID);
		std::vector<esint> partoffset(nextID);
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

		std::vector<std::vector<esint> > frames(nextID), neighbors(nextID);
		#pragma omp parallel for
		for (int p = 0; p < nextID; p++) {
			frames[p].resize(1 + tdecomposition[0][p].size());
			neighbors[p].resize(tdualsize[0][p]);
		}

		// TODO: try parallelization
		// #pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			size_t partindex;
			std::vector<esint> foffset(nextID), noffset(nextID);
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
				for (esint n = frames[partindex][foffset[partindex]] - (dualDist[e + 1] - dualDist[e]); n < frames[partindex][foffset[partindex]]; ++n, ++node) {
					neighbors[partindex][n] = std::lower_bound(tdecomposition[0][partindex].begin(), tdecomposition[0][partindex].end(), *node) - tdecomposition[0][partindex].begin();
				}
			}
		}

		std::vector<esint> pparts(nextID);

		double averageDomainSize = _mesh->elements->size / (double)parts;
		size_t partsCounter = 0;
		for (int p = 0; p < nextID; p++) {
			partsCounter += pparts[p] = std::ceil((frames[p].size() - 1) / averageDomainSize);
			clusters.resize(partsCounter, p);
		}
		_mesh->elements->nclusters = nextID;

		finish("process non-continuous dual graph");

		start("METIS::KWay");
		#pragma omp parallel for
		for (int p = 0; p < nextID; p++) {
			METIS::call(
					run::ecf->decomposition.metis_options,
					frames[p].size() - 1,
					frames[p].data(), neighbors[p].data(),
					0, NULL, NULL,
					pparts[p], partition.data() + partoffset[p]);
		}
		finish("METIS::KWay");

		start("reindex METIS output");

		std::vector<esint> ppartition = partition;
		nextID = 0;
		for (size_t p = 0; p < tdecomposition[0].size(); p++) {
			for (size_t i = 0; i < tdecomposition[0][p].size(); ++i) {
				partition[tdecomposition[0][p][i]] = ppartition[partoffset[p] + i] + nextID;
			}
			nextID += pparts[p];
		}

		finish("reindex METIS output");
	}

	start("post-process domains");

	std::vector<esint> permutation(partition.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (esint i, esint j) {
		if (partition[i] == partition[j]) {
			return i < j;
		}
		return partition[i] < partition[j];
	});

	std::vector<size_t> domainDistribution;
	std::vector<size_t> tdistribution;

	esint partindex = 0;
	auto begin = permutation.begin();
	while (begin != permutation.end()) {
		domainDistribution.push_back(begin - permutation.begin());
		begin = std::lower_bound(begin, permutation.end(), ++partindex, [&] (esint i, esint val) {
			return partition[i] < val;
		});
	}
	domainDistribution.push_back(permutation.size());

	// TODO: improve domain distribution for more complicated decomposition
	if (domainDistribution.size() == threads + 1) {
		tdistribution = std::vector<size_t>(domainDistribution.begin(), domainDistribution.end());
	} else {
		if (domainDistribution.size() < threads + 1) {
			tdistribution = tarray<size_t>::distribute(threads, permutation.size());
		} else {
			double averageThreadSize = _mesh->elements->size / (double)threads;
			tdistribution.push_back(0);
			for (size_t t = 0; t < threads - 1; t++) {
				auto more = std::lower_bound(domainDistribution.begin(), domainDistribution.end(), tdistribution.back() + averageThreadSize);
				if (more == domainDistribution.end()) {
					tdistribution.push_back(domainDistribution.back());
				} else {
					auto less = more - 1;
					if (std::fabs(*less - averageThreadSize * (t + 1)) < std::fabs(*more - averageThreadSize * (t + 1))) {
						tdistribution.push_back(*less);
					} else {
						tdistribution.push_back(*more);
					}
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

	std::vector<esint> domainCounter(threads);
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
	Communication::exscan(_mesh->elements->firstDomain);
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

	finish("post-process domains");

	arrangeElementsPermutation(permutation);
	this->permuteElements(permutation, tdistribution);

	arrangeNodes();
}

void MeshPreprocessing::exchangeElements(const std::vector<esint> &partition)
{
	if (_mesh->nodes->elements == NULL) {
		// need for correctly update nodes ranks
		this->linkNodesAndElements();
	}

	if (_mesh->elements->neighbors == NULL) {
		this->computeElementsNeighbors();
	}

	start("exchanging elements");

	TimeEval timing("EXCHANGE ELEMENTS");
	timing.totalTime.startWithBarrier();

	TimeEvent e1("EE COMPUTE TARGETS"); e1.start();

	// 0. Compute targets
	// 1. Serialize element data
	// 2. Serialize node data
	// 3. Send serialized data to target (new MPI processes)
	// 4. Deserialize data
	// 5. Balance nodes data to threads
	// 6. Re-index elements (IDs have to be always increasing)

	size_t threads = environment->OMP_NUM_THREADS;

	// Step 0: Compute targets

	std::vector<int> targets;
	{ // get target processes
		std::vector<std::vector<int> > ttargets(threads);
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {

			std::vector<int> ttt;
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

	std::vector<std::vector<esint> >  elemsIDs(threads);
	std::vector<std::vector<int> >      elemsBody(threads);
	std::vector<std::vector<int> >      elemsMaterial(threads);
	std::vector<std::vector<Element*> > elemsEpointer(threads);
	std::vector<std::vector<esint> >  elemsNodesDistribution(threads);
	std::vector<std::vector<esint> >  elemsNodesData(threads);
	std::vector<std::vector<esint> >  elemsNeighborsDistribution(threads);
	std::vector<std::vector<esint> >  elemsNeighborsData(threads);
	std::vector<std::vector<esint> >  elemsRegions(threads);

	NodeStore *nodes = new NodeStore();

	std::vector<std::vector<esint> >  nodesIDs(threads);
	std::vector<std::vector<Point> >    nodesCoordinates(threads);
	std::vector<std::vector<esint> >  nodesElemsDistribution(threads);
	std::vector<std::vector<esint> >  nodesElemsData(threads);
	std::vector<std::vector<esint> >  nodesRegions(threads);

	std::vector<std::vector<std::vector<esint> > >  boundaryEDistribution(_mesh->boundaryRegions.size(), std::vector<std::vector<esint> >(threads));
	std::vector<std::vector<std::vector<esint> > >  boundaryEData(_mesh->boundaryRegions.size(), std::vector<std::vector<esint> >(threads));
	std::vector<std::vector<std::vector<Element*> > > boundaryEPointers(_mesh->boundaryRegions.size(), std::vector<std::vector<Element*> >(threads));

	// regions are transfered via mask
	int eregionsBitMaskSize = _mesh->elementsRegions.size() / (8 * sizeof(esint)) + (_mesh->elementsRegions.size() % (8 * sizeof(esint)) ? 1 : 0);
	int bregionsBitMaskSize = _mesh->boundaryRegions.size() / (8 * sizeof(esint)) + (_mesh->boundaryRegions.size() % (8 * sizeof(esint)) ? 1 : 0);

	// serialize data that have to be exchanged
	// the first thread value denotes the thread data size

	// threads x target x elements(id, body, material, code, dualsize, dualdata, nodesize, nodeindices)
	std::vector<std::vector<std::vector<esint> > > sElements(threads, std::vector<std::vector<esint> >(targets.size(), std::vector<esint>({ 0 })));
	std::vector<std::vector<esint> > rElements;

	// threads x target x nodes(id, point, linksize, links, regionMask) + size
	std::vector<std::vector<std::vector<esint> > > sNodes(threads, std::vector<std::vector<esint> >(targets.size(), std::vector<esint>({ 0 })));
	std::vector<std::vector<esint> > rNodes;

	// threads x target x boundary(prefix, (code, nodes))
	std::vector<std::vector<std::vector<esint> > > sBoundary(threads, std::vector<std::vector<esint> >(targets.size(), std::vector<esint>(_mesh->boundaryRegions.size())));
	std::vector<std::vector<esint> > rBoundary;

	e1.end(); timing.addEvent(e1);

	TimeEvent e2("EE SERIALIZE ELEMENTS"); e2.start();

	// Step 1: Serialize element data

	std::vector<esint> regionElementMask(_mesh->elements->size * eregionsBitMaskSize);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		esint maskOffset = 0;
		for (size_t r = 0; r < _mesh->elementsRegions.size(); r++) {
			maskOffset = r / (8 * sizeof(esint));
			esint bit = 1 << (r % (8 * sizeof(esint)));
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
		auto enodes = _mesh->elements->procNodes->cbegin(t);
		auto eneighbors = _mesh->elements->neighbors->cbegin(t);
		auto nIDs = _mesh->nodes->IDs->datatarray().data();

		std::vector<std::vector<esint> > tsElements(targets.size(), std::vector<esint>({ 0 }));

		std::vector<esint>  telemsIDs;
		std::vector<int>      telemsBody;
		std::vector<int>      telemsMaterial;
		std::vector<Element*> telemsEpointer;
		std::vector<esint>  telemsNodesDistribution;
		std::vector<esint>  telemsNodesData;
		std::vector<esint>  telemsNeighborsDistribution;
		std::vector<esint>  telemsNeighborsData;
		std::vector<esint>  telemsRegions;
		if (t == 0) {
			telemsNodesDistribution.push_back(0);
			telemsNeighborsDistribution.push_back(0);
		}

		// estimation
		telemsIDs.reserve(1.5 * _mesh->elements->size / threads);
		telemsBody.reserve(1.5 * _mesh->elements->size / threads);
		telemsMaterial.reserve(1.5 * _mesh->elements->size / threads);
		telemsEpointer.reserve(1.5 * _mesh->elements->size / threads);
		telemsNodesDistribution.reserve(1.5 * _mesh->elements->size / threads);
		telemsNeighborsDistribution.reserve(1.5 * _mesh->elements->size / threads);
		telemsRegions.reserve(1.5 * _mesh->elements->size / threads);

		size_t target;
		for (size_t e = _mesh->elements->distribution[t]; e < _mesh->elements->distribution[t + 1]; ++e, ++enodes, ++eneighbors) {
			if (partition[e] == environment->MPIrank) {
				telemsIDs.push_back(IDs[e]);
				telemsBody.push_back(body[e]);
				telemsMaterial.push_back(material[e]);
				telemsEpointer.push_back(code[e]);
				for (auto n = enodes->begin(); n != enodes->end(); ++n) {
					telemsNodesData.push_back(nIDs[*n]);
				}
				telemsNodesDistribution.push_back(telemsNodesData.size());
				telemsNeighborsData.insert(telemsNeighborsData.end(), eneighbors->begin(), eneighbors->end());
				telemsNeighborsDistribution.push_back(telemsNeighborsData.size());
				telemsRegions.insert(telemsRegions.end(), regionElementMask.begin() + e * eregionsBitMaskSize, regionElementMask.begin() + (e + 1) * eregionsBitMaskSize);
			} else {
				target = t2i(partition[e]);
				tsElements[target].insert(tsElements[target].end(), { IDs[e], body[e], material[e], static_cast<esint>(code[e] - _mesh->_eclasses[t]) });
				tsElements[target].push_back(enodes->size());
				for (auto n = enodes->begin(); n != enodes->end(); ++n) {
					tsElements[target].push_back(nIDs[*n]);
				}
				tsElements[target].push_back(eneighbors->size());
				tsElements[target].insert(tsElements[target].end(), eneighbors->begin(), eneighbors->end());
				tsElements[target].insert(tsElements[target].end(), regionElementMask.begin() + e * eregionsBitMaskSize, regionElementMask.begin() + (e + 1) * eregionsBitMaskSize);
			}
		}

		elemsIDs[t].swap(telemsIDs);
		elemsBody[t].swap(telemsBody);
		elemsMaterial[t].swap(telemsMaterial);
		elemsEpointer[t].swap(telemsEpointer);
		elemsNodesDistribution[t].swap(telemsNodesDistribution);
		elemsNodesData[t].swap(telemsNodesData);
		elemsNeighborsDistribution[t].swap(telemsNeighborsDistribution);
		elemsNeighborsData[t].swap(telemsNeighborsData);
		elemsRegions[t].swap(telemsRegions);

		sElements[t].swap(tsElements);
	}

	e2.end(); timing.addEvent(e2);

	TimeEvent e3("EE SERIALIZE NODES"); e3.start();

	// Step 2: Serialize node data

	std::vector<esint> regionNodeMask(_mesh->nodes->size * bregionsBitMaskSize);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		esint maskOffset = 0;
		for (size_t r = 0; r < _mesh->boundaryRegions.size(); r++) {
			if (_mesh->boundaryRegions[r]->nodes) {
				maskOffset = r / (8 * sizeof(esint));
				esint bit = 1 << (r % (8 * sizeof(esint)));
				auto begin = std::lower_bound(_mesh->boundaryRegions[r]->nodes->datatarray().begin(), _mesh->boundaryRegions[r]->nodes->datatarray().end(), _mesh->nodes->distribution[t]);
				auto end = std::lower_bound(_mesh->boundaryRegions[r]->nodes->datatarray().begin(), _mesh->boundaryRegions[r]->nodes->datatarray().end(), _mesh->nodes->distribution[t + 1]);
				for (auto i = begin; i != end; ++i) {
					regionNodeMask[*i * bregionsBitMaskSize + maskOffset] |= bit;
				}
			}
		}
	}

	esint eBegin = _mesh->elements->gatherElementsProcDistribution()[environment->MPIrank];
	esint eEnd = eBegin + _mesh->elements->size;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		const auto &IDs = _mesh->nodes->IDs->datatarray();
		const auto &coordinates = _mesh->nodes->coordinates->datatarray();
		auto elems = _mesh->nodes->elements->cbegin(t);

		std::vector<esint>  tnodesIDs;
		std::vector<Point>    tnodesCoordinates;
		std::vector<esint>  tnodesElemsDistribution;
		std::vector<esint>  tnodesElemsData;
		std::vector<esint>  tnodesRegions;

		if (t == 0) {
			tnodesElemsDistribution.push_back(0);
		}

		std::vector<std::vector<esint> > tsNodes(targets.size(), std::vector<esint>({ 0 }));

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
						tsNodes[target].insert(tsNodes[target].end(), reinterpret_cast<const esint*>(coordinates.data() + n), reinterpret_cast<const esint*>(coordinates.data() + n + 1));
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

	e3.end(); timing.addEvent(e3);

	TimeEvent e4("EE SERIALIZE BOUNDARY"); e4.start();

	// Step 2.1: Serialize boundary regions data

	std::vector<esint> emembership;

	for (size_t r = 0; r < _mesh->boundaryRegions.size(); r++) {
		if (_mesh->boundaryRegions[r]->dimension) {
			emembership.clear();
			emembership.resize(_mesh->boundaryRegions[r]->distribution.back());
			std::vector<size_t> distribution = _mesh->boundaryRegions[r]->distribution;

			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				auto enodes = _mesh->boundaryRegions[r]->procNodes->cbegin() + distribution[t];
				std::vector<esint> nlinks;
				size_t counter;
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
				std::vector<esint>  tboundaryEDistribution;
				std::vector<esint>  tboundaryEData;
				std::vector<Element*> tboundaryEPointers;
				if (t == 0) {
					tboundaryEDistribution.push_back(0);
				}

				std::vector<size_t> tsize(targets.size());
				for (size_t i = 0; i < targets.size(); i++) {
					tsize[i] = sBoundary[t][i].size();
				}
				esint target;
				auto enodes = _mesh->boundaryRegions[r]->procNodes->cbegin() + distribution[t];
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

	e4.end(); timing.addEvent(e4);

	TimeEvent e5("EE EXCHANGE DATA"); e5.start();

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

	e5.end(); timing.addEvent(e5);

	TimeEvent e6("EE DESERIALIZE ELEMENTS"); e6.start();

	// Step 4: Deserialize element data
	for (size_t i = 0; i < rElements.size(); i++) {
		std::vector<size_t> rdistribution({ 0 });
		size_t p = 0;
		while (p < rElements[i].size()) {
			rdistribution.push_back(p += rElements[i][p]);
		}

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {

			std::vector<esint>  telemsIDs;
			std::vector<int>      telemsBody;
			std::vector<int>      telemsMaterial;
			std::vector<Element*> telemsEpointer;
			std::vector<esint>  telemsNodesDistribution;
			std::vector<esint>  telemsNodesData;
			std::vector<esint>  telemsNeighborsDistribution;
			std::vector<esint>  telemsNeighborsData;
			std::vector<esint>  telemsRegions;

			telemsIDs.reserve(rdistribution[t + 1] - rdistribution[t]);
			telemsBody.reserve(rdistribution[t + 1] - rdistribution[t]);
			telemsMaterial.reserve(rdistribution[t + 1] - rdistribution[t]);
			telemsEpointer.reserve(rdistribution[t + 1] - rdistribution[t]);
			telemsNodesDistribution.reserve(rdistribution[t + 1] - rdistribution[t] + 1);
			telemsNeighborsDistribution.reserve(rdistribution[t + 1] - rdistribution[t] + 1);
			telemsRegions.reserve(rdistribution[t + 1] - rdistribution[t]);

			esint distOffset = 0, neighOffset = 0;

			if (elemsNodesDistribution[t].size()) {
				distOffset = elemsNodesDistribution[t].back();
			}
			if (elemsNeighborsDistribution[t].size()) {
				neighOffset = elemsNeighborsDistribution[t].back();
			}
			if (t == 0 && elemsNodesDistribution[t].size() == 0) {
				telemsNodesDistribution.push_back(0);
			}
			if (t == 0 && elemsNeighborsDistribution[t].size() == 0) {
				telemsNeighborsDistribution.push_back(0);
			}

			for (size_t e = rdistribution[t] + 1; e < rdistribution[t + 1]; ) {
				telemsIDs.push_back(rElements[i][e++]);
				telemsBody.push_back(rElements[i][e++]);
				telemsMaterial.push_back(rElements[i][e++]);
				telemsEpointer.push_back(_mesh->_eclasses[t] + rElements[i][e++]);
				telemsNodesData.insert(telemsNodesData.end(), rElements[i].begin() + e + 1, rElements[i].begin() + e + 1 + rElements[i][e]);
				telemsNodesDistribution.push_back(telemsNodesData.size() + distOffset);
				e += rElements[i][e++]; // nodes + nodes size
				telemsNeighborsData.insert(telemsNeighborsData.end(), rElements[i].begin() + e + 1, rElements[i].begin() + e + 1 + rElements[i][e]);
				telemsNeighborsDistribution.push_back(telemsNeighborsData.size() + neighOffset);
				e += rElements[i][e++]; // neighbors + neighbors size
				telemsRegions.insert(telemsRegions.end(), rElements[i].begin() + e, rElements[i].begin() + e + eregionsBitMaskSize);
				e += eregionsBitMaskSize;
			}

			elemsIDs[t].insert(elemsIDs[t].end(), telemsIDs.begin(), telemsIDs.end());
			elemsBody[t].insert(elemsBody[t].end(), telemsBody.begin(), telemsBody.end());
			elemsMaterial[t].insert(elemsMaterial[t].end(), telemsMaterial.begin(), telemsMaterial.end());
			elemsEpointer[t].insert(elemsEpointer[t].end(), telemsEpointer.begin(), telemsEpointer.end());
			elemsNodesDistribution[t].insert(elemsNodesDistribution[t].end(), telemsNodesDistribution.begin(), telemsNodesDistribution.end());
			elemsNodesData[t].insert(elemsNodesData[t].end(), telemsNodesData.begin(), telemsNodesData.end());
			elemsNeighborsDistribution[t].insert(elemsNeighborsDistribution[t].end(), telemsNeighborsDistribution.begin(), telemsNeighborsDistribution.end());
			elemsNeighborsData[t].insert(elemsNeighborsData[t].end(), telemsNeighborsData.begin(), telemsNeighborsData.end());
			elemsRegions[t].insert(elemsRegions[t].end(), telemsRegions.begin(), telemsRegions.end());
		}
	}

	e6.end(); timing.addEvent(e6);

	TimeEvent e7("EE DESERIALIZE NODES"); e7.start();

	// Step 4: Deserialize node data
	std::vector<esint> nodeset;
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
		std::vector<std::vector<esint> > tnodeset(threads);
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			Point point;

			std::vector<esint>  tnodesIDs;
			std::vector<Point>    tnodesCoordinates;
			std::vector<esint>  tnodesElemsDistribution;
			std::vector<esint>  tnodesElemsData;
			std::vector<esint>  tnodesRegions;
			std::vector<esint>  tnodeSet;
			std::vector<esint>  tnpermutation;

			tnodesIDs.reserve(rdistribution[t + 1] - rdistribution[t]);
			tnodesCoordinates.reserve(rdistribution[t + 1] - rdistribution[t]);
			tnodesElemsDistribution.reserve(rdistribution[t + 1] - rdistribution[t]);
			tnodesRegions.reserve(rdistribution[t + 1] - rdistribution[t]);
			tnodeSet.reserve(rdistribution[t + 1] - rdistribution[t]);

			esint distOffset = 0;

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
				tnpermutation.push_back(n);
				n += 1 + sizeof(Point) / sizeof(esint); // id, Point
				n += 1 + rNodes[i][n]; // linksize, links
				n += bregionsBitMaskSize; // region mask
			}
			std::sort(tnpermutation.begin(), tnpermutation.end(), [&] (esint n1, esint n2) {
				return rNodes[i][n1] < rNodes[i][n2];
			});
			size_t index;
			for (size_t n = 0; n < tnpermutation.size(); n++) {
				index = tnpermutation[n];
				while (nodesetit != nodeset.end() && *nodesetit < rNodes[i][index]) ++nodesetit;
				if (nodesetit == nodeset.end() || *nodesetit != rNodes[i][index]) {
					tnodesIDs.push_back(rNodes[i][index]);
					tnodeSet.push_back(rNodes[i][index]);
					index += 1; //ID
					memcpy(&point, rNodes[i].data() + index, sizeof(Point));
					tnodesCoordinates.push_back(point);
					index += sizeof(Point) / sizeof(esint); // points
					tnodesElemsData.insert(tnodesElemsData.end(), rNodes[i].begin() + index + 1, rNodes[i].begin() + index + 1 + rNodes[i][index]);
					tnodesElemsDistribution.push_back(tnodesElemsData.size() + distOffset);
					index += rNodes[i][index] + 1; // linksize + links
					tnodesRegions.insert(tnodesRegions.end(), rNodes[i].begin() + index, rNodes[i].begin() + index + bregionsBitMaskSize);
					index += bregionsBitMaskSize; // region mask
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
	Esutils::threadDistributionToFullDistribution(elemsNeighborsDistribution);
	Esutils::threadDistributionToFullDistribution(nodesElemsDistribution);

	e7.end(); timing.addEvent(e7);

	TimeEvent e8("EE DESERIALIZE BOUNDARY"); e8.start();

	// Step 4: Deserialize boundary data
	for (size_t n = 0; n < rBoundary.size(); ++n) {
		std::vector<std::vector<esint> > toffset(_mesh->boundaryRegions.size()), tsize(_mesh->boundaryRegions.size());
		esint offset = 0, p = 0;
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

				std::vector<esint>  tboundaryEDistribution;
				std::vector<esint>  tboundaryEData;
				std::vector<Element*> tboundaryEPointers;

				esint distOffset = 0;

				if (boundaryEDistribution[r][t].size()) {
					distOffset = boundaryEDistribution[r][t].back();
				}
				if (t == 0 && boundaryEDistribution[r][t].size() == 0) {
					tboundaryEDistribution.push_back(0);
				}

				for (esint i = toffset[r][t]; i < toffset[r][t] + tsize[r][t];) {
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

	e8.end(); timing.addEvent(e8);

	TimeEvent e9("EE BUILD DATA STORAGES"); e9.start();

	// elements are redistributed later while decomposition -> distribution is not changed now
	std::vector<size_t> elemDistribution(threads + 1);
	for (size_t t = 1; t <= threads; t++) {
		elemDistribution[t] = elemDistribution[t - 1] + elemsIDs[t - 1].size();
	}

	elements->IDs = new serializededata<esint, esint>(1, elemsIDs);
	elements->body = new serializededata<esint, int>(1, elemsBody);
	elements->material = new serializededata<esint, int>(1, elemsMaterial);
	elements->epointers = new serializededata<esint, Element*>(1, elemsEpointer);
	elements->procNodes = new serializededata<esint, esint>(elemsNodesDistribution, elemsNodesData); // global IDs

	elements->regionMaskSize = eregionsBitMaskSize;
	elements->regions = new serializededata<esint, esint>(eregionsBitMaskSize, elemsRegions);

	elements->neighbors = new serializededata<esint, esint>(elemsNeighborsDistribution, elemsNeighborsData);

	elements->dimension = _mesh->elements->dimension;
	elements->size = elements->IDs->structures();
	elements->distribution = elements->IDs->datatarray().distribution();

	// Step 5: Balance node data to threads
	std::vector<size_t> nodeDistribution(threads);
	for (size_t t = 1; t < threads; t++) {
		nodeDistribution[t] = nodeDistribution[t - 1] + nodesIDs[t - 1].size();
	}

	serializededata<esint, esint>::balance(1, nodesIDs);
	serializededata<esint, Point>::balance(1, nodesCoordinates);
	serializededata<esint, esint>::balance(nodesElemsDistribution, nodesElemsData);

	nodes->IDs = new serializededata<esint, esint>(1, nodesIDs);
	nodes->coordinates = new serializededata<esint, Point>(1, nodesCoordinates);
	nodes->elements = new serializededata<esint, esint>(nodesElemsDistribution, nodesElemsData);
	nodes->size = nodes->IDs->datatarray().size();
	nodes->distribution = nodes->IDs->datatarray().distribution();

	for (size_t r = 0; r < _mesh->boundaryRegions.size(); r++) {
		if (_mesh->boundaryRegions[r]->dimension) {
			delete _mesh->boundaryRegions[r]->procNodes;
			delete _mesh->boundaryRegions[r]->epointers;

			Esutils::threadDistributionToFullDistribution(boundaryEDistribution[r]);
			_mesh->boundaryRegions[r]->procNodes = new serializededata<esint, esint>(boundaryEDistribution[r], boundaryEData[r]);
			_mesh->boundaryRegions[r]->epointers = new serializededata<esint, Element*>(1, boundaryEPointers[r]);
		}
	}

	for (size_t r = 0; r < _mesh->elementsRegions.size(); r++) {
		esint maskOffset = r / (8 * sizeof(esint));
		esint bit = 1 << (r % (8 * sizeof(esint)));
		delete _mesh->elementsRegions[r]->elements;
		std::vector<std::vector<esint> > regionelems(threads);

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (size_t i = 0; i < elemsRegions[t].size(); i += eregionsBitMaskSize) {
				if (elemsRegions[t][i + maskOffset] & bit) {
					regionelems[t].push_back(elemDistribution[t] + i / eregionsBitMaskSize);
				}
			}
		}

		serializededata<esint, esint>::balance(1, regionelems);
		_mesh->elementsRegions[r]->elements = new serializededata<esint, esint>(1, regionelems);
	}

	for (size_t r = 0; r < _mesh->boundaryRegions.size(); r++) {
		if (_mesh->boundaryRegions[r]->nodes) {
			esint maskOffset = r / (8 * sizeof(esint));
			esint bit = 1 << (r % (8 * sizeof(esint)));
			delete _mesh->boundaryRegions[r]->nodes;
			std::vector<std::vector<esint> > regionnodes(threads);

			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				for (size_t i = 0; i < nodesRegions[t].size(); i += bregionsBitMaskSize) {
					if (nodesRegions[t][i + maskOffset] & bit) {
						regionnodes[t].push_back(nodeDistribution[t] + i / bregionsBitMaskSize);
					}
				}
			}

			serializededata<esint, esint>::balance(1, regionnodes);
			_mesh->boundaryRegions[r]->nodes = new serializededata<esint, esint>(1, regionnodes);
		}
	}

	e9.end(); timing.addEvent(e9);

	TimeEvent e10("EE SYNCHRONIZE ELEMENTS IDS"); e10.start();

	std::vector<esint> eIDsOLD = _mesh->elements->gatherElementsProcDistribution();
	std::vector<esint> eIDsNEW = elements->gatherElementsProcDistribution();

	for (size_t t = 1; t < threads; ++t) {
		elemsIDs[0].insert(elemsIDs[0].end(), elemsIDs[t].begin(), elemsIDs[t].end());
	}

	std::vector<esint> epermutation(elemsIDs[0].size());
	std::iota(epermutation.begin(), epermutation.end(), 0);
	std::sort(epermutation.begin(), epermutation.end(), [&] (esint i, esint j) { return elemsIDs[0][i] < elemsIDs[0][j]; });

	std::vector<esint> sortedElements(epermutation.size());
	#pragma omp parallel for
	for (size_t t = 0; t < threads; ++t) {
		for (size_t e = elemDistribution[t]; e < elemDistribution[t + 1]; e++) {
			sortedElements[e] = elemsIDs[0][epermutation[e]];
		}
		Esutils::sortAndRemoveDuplicity(nodesElemsData[t]);
	}
	Esutils::inplaceMerge(nodesElemsData);
	Esutils::removeDuplicity(nodesElemsData[0]);

	std::vector<std::vector<esint> > requestedIDs, receivedTargets, IDrequests;
	std::vector<std::vector<int> > IDtargets(threads);
	std::vector<int> sources;

	std::vector<size_t> rdistribution = tarray<size_t>::distribute(threads, environment->MPIsize);

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
		auto mybegin = std::lower_bound(sortedElements.begin(), sortedElements.end(), eIDsOLD[IDtargets[0][t]]);
		auto myend   = std::lower_bound(sortedElements.begin(), sortedElements.end(), eIDsOLD[IDtargets[0][t] + 1]);
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

	std::vector<std::vector<bool> > newtargets(threads, std::vector<bool>(environment->MPIsize, false));
	for (size_t r = 0; r < receivedTargets.size(); r++) {
		rdistribution = tarray<size_t>::distribute(threads, receivedTargets[r].size());

		#pragma omp parallel for
		for (size_t t = 0; t < threads; ++t) {
			for (size_t e = rdistribution[t]; e < rdistribution[t + 1]; ++e) {
				newtargets[t][receivedTargets[r][e]] = true;
			}
		}
	}

	#pragma omp parallel for
	for (size_t t = 0; t < threads; ++t) {
		std::vector<int> ttargets;
		for (int r = 0; r < environment->MPIsize; r++) {
			if (r != environment->MPIrank && newtargets[t][r]) {
				ttargets.push_back(r);
			}
		}
		IDtargets[t].swap(ttargets);
	}
	for (size_t t = 1; t < threads; ++t) {
		IDtargets[0].insert(IDtargets[0].end(), IDtargets[t].begin(), IDtargets[t].end());
	}
	Esutils::sortAndRemoveDuplicity(IDtargets[0]);

	std::vector<std::vector<std::vector<esint> > > newIDrequests(threads, std::vector<std::vector<esint> >(IDtargets[0].size()));
	std::vector<std::vector<esint> > newIDs;

	for (size_t r = 0; r < receivedTargets.size(); r++) {
		rdistribution = tarray<size_t>::distribute(threads, receivedTargets[r].size());

		#pragma omp parallel for
		for (size_t t = 0; t < threads; ++t) {
			size_t tindex;
			std::vector<std::vector<esint> > tnewIDrequests(IDtargets[0].size());

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
				auto eit = std::lower_bound(sortedElements.begin(), sortedElements.end(), IDrequests[r][rdistribution[t]]);
				for (size_t e = rdistribution[t]; e < rdistribution[t + 1]; ++e) {
					while (eit != sortedElements.end() && *eit < IDrequests[r][e]) ++eit;
					IDrequests[r][e] = epermutation[eit - sortedElements.begin()] + eIDsNEW[environment->MPIrank];
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
	std::vector<std::pair<esint, esint> > IDMap(offsets.back());

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

	e10.end(); timing.addEvent(e10);

	TimeEvent e11("EE REINDEX"); e11.start();

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (auto elem = nodes->elements->begin(t); elem != nodes->elements->end(t); ++elem) {
			for (auto e = elem->begin(); e != elem->end(); ++e) {
				auto mapit = std::lower_bound(IDMap.begin(), IDMap.end(), std::make_pair(*e, (esint)0));
				if (mapit == IDMap.end() || mapit->first != *e) {
					*e = epermutation[std::lower_bound(sortedElements.begin(), sortedElements.end(), *e) - sortedElements.begin()] + eIDsNEW[environment->MPIrank];
				} else {
					*e = mapit->second;
				}
			}
			std::sort(elem->begin(), elem->end());
		}
	}

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (auto neighbors = elements->neighbors->begin(t); neighbors != elements->neighbors->end(t); ++neighbors) {
			for (auto n = neighbors->begin(); n != neighbors->end(); ++n) {
				if (*n != -1) {
					auto mapit = std::lower_bound(IDMap.begin(), IDMap.end(), std::make_pair(*n, (esint)0));
					if (mapit == IDMap.end() || mapit->first != *n) {
						*n = epermutation[std::lower_bound(sortedElements.begin(), sortedElements.end(), *n) - sortedElements.begin()] + eIDsNEW[environment->MPIrank];
					} else {
						*n = mapit->second;
					}
				}
			}
		}
	}


	std::vector<esint> permutation(nodes->size);
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (esint i, esint j) { return nodes->IDs->datatarray()[i] < nodes->IDs->datatarray()[j]; });

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (auto n = elements->procNodes->begin(t)->begin(); n != elements->procNodes->end(t)->begin(); ++n) {
			*n = *std::lower_bound(permutation.begin(), permutation.end(), *n, [&] (esint i, esint val) {
				return nodes->IDs->datatarray()[i] < val;
			});
		}
	}

	for (size_t r = 0; r < _mesh->boundaryRegions.size(); r++) {
		if (_mesh->boundaryRegions[r]->procNodes) {
			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				for (auto n = _mesh->boundaryRegions[r]->procNodes->begin(t)->begin(); n != _mesh->boundaryRegions[r]->procNodes->end(t)->begin(); ++n) {
					*n = *std::lower_bound(permutation.begin(), permutation.end(), *n, [&] (esint i, esint val) {
						return nodes->IDs->datatarray()[i] < val;
					});
				}
			}
		}
	}

	std::vector<std::vector<esint> > rankBoundaries(threads);
	std::vector<std::vector<int> > rankData(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<esint> trankBoundaries;
		std::vector<int> trankData;
		if (t == 0) {
			trankBoundaries.push_back(0);
		}

		int rank;
		for (auto elem = nodes->elements->begin(t); elem != nodes->elements->end(t); ++elem) {
			trankData.push_back(std::lower_bound(eIDsNEW.begin(), eIDsNEW.end(), *elem->begin() + 1) - eIDsNEW.begin() - 1);
			for (auto e = elem->begin() + 1; e != elem->end(); ++e) {
				rank = std::lower_bound(eIDsNEW.begin(), eIDsNEW.end(), *e + 1) - eIDsNEW.begin() - 1;
				if (rank != trankData.back()) {
					trankData.push_back(rank);
				}
			}
			trankBoundaries.push_back(trankData.size());
		}

		rankBoundaries[t].swap(trankBoundaries);
		rankData[t].swap(trankData);
	}

	Esutils::threadDistributionToFullDistribution(rankBoundaries);

	nodes->ranks = new serializededata<esint, int>(rankBoundaries, rankData);

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

	e11.end(); timing.addEvent(e11);

	timing.totalTime.endWithBarrier();
	timing.printStatsMPI();

	delete elements;
	delete nodes;

	finish("exchanging elements");
}

void MeshPreprocessing::permuteElements(const std::vector<esint> &permutation, const std::vector<size_t> &distribution)
{
	if (_mesh->nodes->elements == NULL) {
		this->linkNodesAndElements();
	}

	start("permutation of elements");

	std::vector<esint> backpermutation(permutation.size());
	std::iota(backpermutation.begin(), backpermutation.end(), 0);
	std::sort(backpermutation.begin(), backpermutation.end(), [&] (esint i, esint j) { return permutation[i] < permutation[j]; });

	size_t threads = environment->OMP_NUM_THREADS;

	auto n2i = [ & ] (size_t neighbor) {
		return std::lower_bound(_mesh->neighbours.begin(), _mesh->neighbours.end(), neighbor) - _mesh->neighbours.begin();
	};

	std::vector<esint> IDBoundaries = _mesh->elements->gatherElementsProcDistribution();
	std::vector<std::vector<std::pair<esint, esint> > > rHalo(_mesh->neighbours.size());

	if (_mesh->elements->neighbors != NULL || _mesh->nodes->elements != NULL) {
		// thread x neighbor x elements(oldID, newID)
		std::vector<std::vector<std::vector<std::pair<esint, esint> > > > sHalo(threads, std::vector<std::vector<std::pair<esint, esint> > >(_mesh->neighbours.size()));

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			auto ranks = _mesh->nodes->ranks->cbegin(t);
			auto elements = _mesh->nodes->elements->cbegin(t);
			esint begine = IDBoundaries[environment->MPIrank];
			esint ende   = IDBoundaries[environment->MPIrank + 1];

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

	auto globalremap = [&] (serializededata<esint, esint>* data, bool sort) {
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
							*n = std::lower_bound(rHalo[n2i(source)].begin(), rHalo[n2i(source)].end(), std::make_pair(*n, (esint)0))->second;
						}
					}
				}
				if (sort) {
					std::sort(e->begin(), e->end());
				}
			}
		}
	};

	esint firstID = _mesh->elements->IDs->datatarray().front();
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
