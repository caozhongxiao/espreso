
#include "transformations.h"

#include "../newmesh.h"
#include "../elements/elementstore.h"

#include "../../basis/containers/serializededata.h"

#include "../../basis/logging/logging.h"
#include "../../basis/utilities/utils.h"
#include "../../basis/utilities/communication.h"

#include "../../wrappers/wmetis.h"
#include "../../wrappers/wparmetis.h"

#include "../../config/ecf/environment.h"

using namespace espreso;


void Transformation::reclusterize(NewMesh &mesh)
{
	ESINFO(TVERBOSITY) << std::string(2 * level++, ' ') << "Transformation::re-distribution of the mesh to processes started.";

	if (mesh._elems->dual == NULL) {
		Transformation::computeDual(mesh);
	}
	if (mesh._elems->coordinates == NULL) {
		Transformation::computeElementCenters(mesh);
	}

	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<esglobal> edistribution(environment->MPIsize + 1);

	esglobal esize = mesh._elems->size;
	Communication::exscan(esize);

	MPI_Allgather(&esize, sizeof(esglobal), MPI_BYTE, edistribution.data(), sizeof(esglobal), MPI_BYTE, MPI_COMM_WORLD);
	edistribution.back() = esize + mesh._elems->size;
	MPI_Bcast(&edistribution.back(), sizeof(esglobal), MPI_BYTE, environment->MPIsize - 1, MPI_COMM_WORLD);

	if (environment->MPIsize == 1) {
		ESINFO(TVERBOSITY) << "Transformation::re-distribution of the mesh to processes skipped (there is only 1 MPI process).";
		return;
	}


	std::vector<esglobal> partition(mesh._elems->size), permutation(mesh._elems->size), edgeWeights(mesh._elems->dual->data().size());

	size_t edgeConst = 10000;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto dual = mesh._elems->dual->cbegin(t);
		int material;
		NewElement::TYPE type;

		for (size_t e = mesh._elems->distribution[t]; e < mesh._elems->distribution[t + 1]; ++e, ++dual) {
			for (auto neigh = dual->begin(); neigh != dual->end(); ++neigh) {
				auto it = std::lower_bound(mesh._elems->IDs->data().cbegin(), mesh._elems->IDs->data().cend(), *neigh);
				if (it == mesh._elems->IDs->data().cend()) {
					auto halo = std::lower_bound(mesh._halo->IDs->data().cbegin(), mesh._halo->IDs->data().cend(), *neigh);
					material = mesh._halo->material->data()[halo - mesh._halo->IDs->data().cbegin()];
					type = mesh._halo->epointers->data()[halo - mesh._halo->IDs->data().cbegin()]->type;
				} else {
					material = mesh._elems->material->data()[it - mesh._elems->IDs->data().cbegin()];
					type = mesh._elems->epointers->data()[it - mesh._elems->IDs->data().cbegin()]->type;
				}
			}
		}

//		while (dual.next()) {
//			body1 = (*mesh.__elements->elements)[dual.index()]->param(Element::Params::BODY);
//			etype1 = (int)(*mesh.__elements->elements)[dual.index()]->type();
//			material1 = (*mesh.__elements->elements)[dual.index()]->param(Element::Params::MATERIAL);
//			for (auto n = dual.begin(); n != dual.end(); ++n) {
//				if (mesh.__elements->elementOffset <= *n && *n < mesh.__elements->elementOffset + mesh.__elements->elements->size()) {
//					body2 = (*mesh.__elements->elements)[*n - mesh.__elements->elementOffset]->param(Element::Params::BODY);
//					etype2 = (int)(*mesh.__elements->elements)[*n - mesh.__elements->elementOffset]->type();
//					material2 = (*mesh.__elements->elements)[*n - mesh.__elements->elementOffset]->param(Element::Params::MATERIAL);
//				} else {
//					tmph.id = *n;
//					he = &*std::lower_bound(mesh.__elements->haloElements->begin(), mesh.__elements->haloElements->end(), tmph);
//					body2 = he->body;
//					etype2 = he->type;
//					material2 = he->material;
//				}
//
//				edgeWeights[edgeIndex] = 6 * edgeConst + 1;
//				if (body1 != body2) {
//					edgeWeights[edgeIndex] -= 3 * edgeConst;
//				}
//				if (etype1 != etype2) {
//					edgeWeights[edgeIndex] -= 2 * edgeConst;
//				}
//				if (material1 != material2) {
//					edgeWeights[edgeIndex] -= 1 * edgeConst;
//				}
//				edgeIndex++;
//			}
//		}
	}

//	esglobal edgecut = ParMETIS::call(ParMETIS::METHOD::ParMETIS_V3_PartKway,
//		edistribution.data(),
//		mesh.__elements->fullDual->crs(), mesh.__elements->fullDual->data(),
//		mesh.__elements->dimension, mesh.__elements->coordinates->data(),
//		0, NULL, edgeWeights,
//		partition
//	);
//
//	ESINFO(TVERBOSITY) << Info::plain() << "Using ParMETIS to improve edge-cuts: " << edgecut;
//	esglobal prev = 2 * edgecut;
//	while (1.01 * edgecut < prev) {
//		prev = edgecut;
//		edgecut = ParMETIS::call(ParMETIS::METHOD::ParMETIS_V3_AdaptiveRepart,
//			edistribution.data(),
//			mesh.__elements->fullDual->crs(), mesh.__elements->fullDual->data(),
//			mesh.__elements->dimension, mesh.__elements->coordinates->data(),
//			0, NULL, edgeWeights,
//			partition
//		);
//		ESINFO(TVERBOSITY) << Info::plain() << " -> " << edgecut;
//	}
//	ESINFO(TVERBOSITY);
//	delete[] edgeWeights;
//
//	ElementStore *repartition = new ElementStore();
//
//	repartition->distribution  = mesh.__elements->distribution;
//	repartition->dimension     = mesh.__elements->dimension;
//	repartition->neighbors     = mesh.neighbours();
//	repartition->elementOffset = mesh.__elements->elementOffset;
//	repartition->fullDual      = mesh.__elements->fullDual;
//
//	_distributeDualGraph(mesh, edistribution, repartition, partition, permutation);
////	while (!_checkContinuity(mesh)) {
////		_tryrepartition(mesh, permutation);
////		_distributeDualGraph(mesh, partition, permutation);
////	}
//
//	_distributeNewMesh(mesh, repartition, permutation);
//	for (size_t e = 0; e < mesh.__elements->elements->size(); e++) {
//		delete (*mesh.__elements->elements)[e];
//	}
//	delete mesh.__elements;
//	mesh.__elements = repartition;
//
//	delete[] partition;
//	delete[] permutation;

	ESINFO(TVERBOSITY) << std::string(--level * 2, ' ') << "Transformation::re-distribution of the mesh to processes finished.";

//	if (mesh.__elements->fullDual == NULL) {
//		Transformation::createFullDualGraph(mesh);
//	}
//	if (mesh.__elements->coordinates == NULL) {
//		Transformation::computeElementCenters(mesh);
//	}
//
//	std::vector<esglobal> edistribution(environment->MPIsize + 1);
//	MPI_Allgather(&mesh.__elements->elementOffset, sizeof(esglobal), MPI_BYTE, edistribution.data(), sizeof(esglobal), MPI_BYTE, MPI_COMM_WORLD);
//	edistribution.back() = mesh.__elements->elementOffset + mesh.__elements->elements->size();
//	MPI_Bcast(&edistribution.back(), sizeof(esglobal), MPI_BYTE, environment->MPIsize - 1, MPI_COMM_WORLD);
//
//	if (environment->MPIsize == 1) {
//		if (mesh.__elements->intervals.size() < 1) {
//			mesh.__elements->intervals.push_back(std::vector<size_t>(edistribution.begin(), edistribution.end()));
//		} else {
//			mesh.__elements->intervals[0] = std::vector<size_t>(edistribution.begin(), edistribution.end());
//		}
//		ESINFO(TVERBOSITY) << "Transformation::mesh distribution to MPI processes skipped.";
//		return;
//	}
//
//	ESINFO(TVERBOSITY) << "Transformation::mesh distribution to MPI processes started.";
//
//	esglobal *partition   = new esglobal[mesh.__elements->elements->size()];
//	esglobal *permutation = new esglobal[mesh.__elements->elements->size()];
//	esglobal *edgeWeights = new esglobal[mesh.__elements->fullDual->pdata().size()];
//
//	size_t edgeConst = 10000;
//
//	#pragma omp parallel for
//	for (size_t t = 0; t < mesh.__elements->fullDual->pdata().threads(); t++) {
//		const crsiterator<esglobal> &dual = mesh.__elements->fullDual->iterator(t);
//		size_t edgeIndex = mesh.__elements->fullDual->pdata().distribution()[t];
//		int material1, material2, body1, body2, etype1, etype2;
//		HaloElement tmph, *he;
//		while (dual.next()) {
//			body1 = (*mesh.__elements->elements)[dual.index()]->param(Element::Params::BODY);
//			etype1 = (int)(*mesh.__elements->elements)[dual.index()]->type();
//			material1 = (*mesh.__elements->elements)[dual.index()]->param(Element::Params::MATERIAL);
//			for (auto n = dual.begin(); n != dual.end(); ++n) {
//				if (mesh.__elements->elementOffset <= *n && *n < mesh.__elements->elementOffset + mesh.__elements->elements->size()) {
//					body2 = (*mesh.__elements->elements)[*n - mesh.__elements->elementOffset]->param(Element::Params::BODY);
//					etype2 = (int)(*mesh.__elements->elements)[*n - mesh.__elements->elementOffset]->type();
//					material2 = (*mesh.__elements->elements)[*n - mesh.__elements->elementOffset]->param(Element::Params::MATERIAL);
//				} else {
//					tmph.id = *n;
//					he = &*std::lower_bound(mesh.__elements->haloElements->begin(), mesh.__elements->haloElements->end(), tmph);
//					body2 = he->body;
//					etype2 = he->type;
//					material2 = he->material;
//				}
//
//				edgeWeights[edgeIndex] = 6 * edgeConst + 1;
//				if (body1 != body2) {
//					edgeWeights[edgeIndex] -= 3 * edgeConst;
//				}
//				if (etype1 != etype2) {
//					edgeWeights[edgeIndex] -= 2 * edgeConst;
//				}
//				if (material1 != material2) {
//					edgeWeights[edgeIndex] -= 1 * edgeConst;
//				}
//				edgeIndex++;
//			}
//		}
//	}
//
//	esglobal edgecut = ParMETIS::call(ParMETIS::METHOD::ParMETIS_V3_PartKway,
//		edistribution.data(),
//		mesh.__elements->fullDual->crs(), mesh.__elements->fullDual->data(),
//		mesh.__elements->dimension, mesh.__elements->coordinates->data(),
//		0, NULL, edgeWeights,
//		partition
//	);
//
//	ESINFO(TVERBOSITY) << Info::plain() << "Using ParMETIS to improve edge-cuts: " << edgecut;
//	esglobal prev = 2 * edgecut;
//	while (1.01 * edgecut < prev) {
//		prev = edgecut;
//		edgecut = ParMETIS::call(ParMETIS::METHOD::ParMETIS_V3_AdaptiveRepart,
//			edistribution.data(),
//			mesh.__elements->fullDual->crs(), mesh.__elements->fullDual->data(),
//			mesh.__elements->dimension, mesh.__elements->coordinates->data(),
//			0, NULL, edgeWeights,
//			partition
//		);
//		ESINFO(TVERBOSITY) << Info::plain() << " -> " << edgecut;
//	}
//	ESINFO(TVERBOSITY);
//	delete[] edgeWeights;
//
//	ElementStore *repartition = new ElementStore();
//
//	repartition->distribution  = mesh.__elements->distribution;
//	repartition->dimension     = mesh.__elements->dimension;
//	repartition->neighbors     = mesh.neighbours();
//	repartition->elementOffset = mesh.__elements->elementOffset;
//	repartition->fullDual      = mesh.__elements->fullDual;
//
//	_distributeDualGraph(mesh, edistribution, repartition, partition, permutation);
////	while (!_checkContinuity(mesh)) {
////		_tryrepartition(mesh, permutation);
////		_distributeDualGraph(mesh, partition, permutation);
////	}
//
//	_distributeNewMesh(mesh, repartition, permutation);
//	for (size_t e = 0; e < mesh.__elements->elements->size(); e++) {
//		delete (*mesh.__elements->elements)[e];
//	}
//	delete mesh.__elements;
//	mesh.__elements = repartition;
//
//	delete[] partition;
//	delete[] permutation;
//
//	ESINFO(TVERBOSITY) << "Transformation::mesh distribution to MPI processes finished.";
}

void Transformation::partitiate(NewMesh &mesh, esglobal parts, TFlags::SEPARATE separate)
{
//	if (mesh.__elements->restrictedDual == NULL) {
//		Transformation::createRestrictedDualGraph(mesh, separate);
//	}
//
//	if (parts == 1) {
//		if (mesh.__elements->intervals.size() <= mesh.__elements->levelFETI[1]) {
//			mesh.__elements->intervals.resize(mesh.__elements->levelFETI[1] + 1);
//		}
//		mesh.__elements->intervals[mesh.__elements->levelFETI[1]] = mesh.__elements->intervals[mesh.__elements->levelFETI[0]];
//		ESINFO(TVERBOSITY) << "Transformation::mesh partition to domains skipped.";
//		return;
//	}
//	if (parts < 1) {
//		ESINFO(ERROR) << "Invalid number of parts requested.";
//	}
//	ESINFO(TVERBOSITY) << "Transformation::mesh partition to domains started.";
//
//	size_t threads = environment->OMP_NUM_THREADS;
//
//	size_t esize = mesh.__elements->elements->size();
//	size_t eoffset = mesh.__elements->elementOffset;
//	size_t averageSize = esize / parts;
//	size_t edgeConst = 10000;
//
//	esglobal *partition   = new esglobal[esize];
//	esglobal *permutation = new esglobal[esize];
//
//	std::vector<bool> processed(esize, false);
//	size_t partsSum = 0;
//
//	std::vector<eslocal> chunk;
//	std::vector<eslocal> stack;
//	const crsiterator<eslocal> dual = mesh.__elements->restrictedDual->iterator();
//
//	for (size_t i = 0; i < processed.size(); i++) {
//		if (!processed[i]) {
//
//			chunk.clear();
//			stack.clear();
//			stack.push_back(i);
//			processed[i] = true;
//
//			while (stack.size()) {
//				chunk.push_back(stack.back());
//				dual.seek(stack.back());
//				stack.pop_back();
//				for (auto n = dual.begin(); n != dual.end(); ++n) {
//					if (!processed[*n]) {
//						stack.push_back(*n);
//						processed[*n] = true;
//					}
//				}
//			}
//
//			std::sort(chunk.begin(), chunk.end());
//
//			if (chunk.size() == processed.size()) {
//				// elements are continuous
//
//				esglobal *edgeWeights = new esglobal[mesh.__elements->restrictedDual->pdata().size()];
//
//				#pragma omp parallel for
//				for (size_t t = 0; t < mesh.__elements->restrictedDual->pdata().threads(); t++) {
//					const crsiterator<esglobal> &tdual = mesh.__elements->restrictedDual->iterator(t);
//					size_t edgeIndex = mesh.__elements->restrictedDual->pdata().distribution()[t];
//					int material1, material2, body1, body2, etype1, etype2;
//					while (tdual.next()) {
//						body1 = (*mesh.__elements->elements)[tdual.index()]->param(Element::Params::BODY);
//						etype1 = (int)(*mesh.__elements->elements)[tdual.index()]->type();
//						material1 = (*mesh.__elements->elements)[tdual.index()]->param(Element::Params::MATERIAL);
//						for (auto n = tdual.begin(); n != tdual.end(); ++n) {
//							body2 = (*mesh.__elements->elements)[*n]->param(Element::Params::BODY);
//							etype2 = (int)(*mesh.__elements->elements)[*n]->type();
//							material2 = (*mesh.__elements->elements)[*n]->param(Element::Params::MATERIAL);
//
//							edgeWeights[edgeIndex] = 6 * edgeConst + 1;
//							if (body1 != body2) {
//								edgeWeights[edgeIndex] -= 3 * edgeConst;
//							}
//							if (etype1 != etype2) {
//								edgeWeights[edgeIndex] -= 2 * edgeConst;
//							}
//							if (material1 != material2) {
//								edgeWeights[edgeIndex] -= 1 * edgeConst;
//							}
//							edgeIndex++;
//						}
//					}
//				}
//
//				METIS::call(chunk.size(), mesh.__elements->restrictedDual->crs(), mesh.__elements->restrictedDual->data(), 0, NULL, edgeWeights, parts, partition);
//
//				partsSum += parts;
//				delete[] edgeWeights;
//				break;
//			} else {
//				std::sort(chunk.begin(), chunk.end());
//				size_t chunkParts = 1 + (chunk.size() - 1) / averageSize;
//
//				std::vector<size_t> distribution = Esutils::getDistribution(threads, chunk.size());
//				std::vector<size_t> offsets(threads);
//				#pragma omp parallel for
//				for (size_t t = 0; t < threads; t++) {
//					size_t offset = 0;
//					const crsiterator<esglobal> tdual = mesh.__elements->restrictedDual->iterator();
//					for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {
//						tdual.seek(chunk[e]);
//						offset += tdual.size();
//					}
//					offsets[t] = offset;
//				}
//				size_t size = Esutils::sizesToOffsets(offsets);
//				esglobal *crs = new esglobal[chunk.size() + 1];
//				esglobal *data = new esglobal[size];
//				esglobal *edgeWeights = new esglobal[size];
//
//				crs[0] = 0;
//				#pragma omp parallel for
//				for (size_t t = 0; t < threads; t++) {
//					size_t offset = offsets[t];
//					int material1, material2, body1, body2, etype1, etype2;
//					const crsiterator<esglobal> tdual = mesh.__elements->restrictedDual->iterator();
//					for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {
//						tdual.seek(chunk[e]);
//
//						body1 = (*mesh.__elements->elements)[chunk[e]]->param(Element::Params::BODY);
//						etype1 = (int)(*mesh.__elements->elements)[chunk[e]]->type();
//						material1 = (*mesh.__elements->elements)[chunk[e]]->param(Element::Params::MATERIAL);
//						for (auto n = tdual.begin(); n != tdual.end(); ++n) {
//							body2 = (*mesh.__elements->elements)[*n]->param(Element::Params::BODY);
//							etype2 = (int)(*mesh.__elements->elements)[*n]->type();
//							material2 = (*mesh.__elements->elements)[*n]->param(Element::Params::MATERIAL);
//
//							edgeWeights[offset] = 6 * edgeConst + 1;
//							if (body1 != body2) {
//								edgeWeights[offset] -= 3 * edgeConst;
//							}
//							if (etype1 != etype2) {
//								edgeWeights[offset] -= 2 * edgeConst;
//							}
//							if (material1 != material2) {
//								edgeWeights[offset] -= 1 * edgeConst;
//							}
//							data[offset++] = std::lower_bound(chunk.begin(), chunk.end(), *n) - chunk.begin();
//						}
//						crs[e + 1] = offset;
//					}
//				}
//
//				esglobal *subPartition = new esglobal[chunk.size()];
//				METIS::call(chunk.size(), crs, data, 0, NULL, edgeWeights, chunkParts, subPartition);
//				for (size_t e = 0; e < chunk.size(); e++) {
//					partition[chunk[e]] = subPartition[e] + partsSum;
//				}
//
//				partsSum += chunkParts;
//				delete[] crs;
//				delete[] data;
//				delete[] edgeWeights;
//				delete[] subPartition;
//			}
//		}
//	}
//
//	std::iota(permutation, permutation + esize, 0);
//	std::sort(permutation, permutation + esize, [&] (esglobal i, esglobal j) { return partition[i] < partition[j]; });
//	Transformation::_permuteNewMesh(mesh, permutation);
//
//	std::vector<size_t> intervals = { eoffset };
//	for (size_t p = 0; p < partition[permutation[esize - 1]]; p++) {
//		intervals.push_back(eoffset + std::lower_bound(permutation, permutation + esize, p + 1, [&] (esglobal i, esglobal part) { return partition[i] < part; }) - permutation);
//	}
//	if (environment->MPIrank + 1 == environment->MPIsize) {
//		intervals.push_back(eoffset + esize);
//	}
//	Esutils::removeDuplicity(intervals);
//
//	if (mesh.__elements->intervals.size() <= mesh.__elements->levelFETI[1]) {
//		mesh.__elements->intervals.resize(mesh.__elements->levelFETI[1] + 1);
//	}
//	if (!Communication::allgatherUnknownSize(intervals, mesh.__elements->intervals[mesh.__elements->levelFETI[1]])) {
//		ESINFO(ERROR) << "ESPRESO internal error while gathering domain intervals.";
//	}
//
//	delete[] partition;
//	delete[] permutation;
//	ESINFO(TVERBOSITY) << "Transformation::mesh partition to domains finished.";
}


void Transformation::_distributeDualGraph(NewMesh &mesh, std::vector<esglobal> &edistribution, ElementStore *store, esglobal *partition, esglobal *permutation)
{
//	ESINFO(TVERBOSITY) << "Transformation::  distribution of Dual Graph started.";
//
//	size_t threads = environment->OMP_NUM_THREADS;
//	// thread x MPIrank
//	std::vector<std::vector<esglobal> > tsize(threads, std::vector<esglobal>(environment->MPIsize));
//	std::vector<esglobal> csize(environment->MPIsize);
//	std::vector<std::vector<int> > ttargets(threads);
//	std::vector<int> targets;
//
//	#pragma omp parallel for
//	for (size_t t = 0; t < threads; t++) {
//		for (size_t e = store->distribution[t]; e < store->distribution[t + 1]; e++) {
//			tsize[t][partition[e]]++;
//		}
//	}
//
//	std::vector<int> rdistribution = Esutils::getDistribution(threads, environment->MPIsize);
//
//	#pragma omp parallel for
//	for (size_t t = 0; t < threads; t++) {
//		for (int r = rdistribution[t]; r < rdistribution[t + 1]; r++) {
//			for (size_t i = 0; i < threads; i++) {
//				csize[r] += tsize[i][r];
//			}
//			if (csize[r]) {
//				ttargets[t].push_back(r);
//			}
//		}
//	}
//
//	for (size_t t = 0; t < threads; t++) {
//		targets.insert(targets.end(), ttargets[t].begin(), ttargets[t].end());
//	}
//
//	std::vector<esglobal> coffset = Communication::exscan(csize);
//	std::vector<esglobal> roffset(coffset);
//
//	std::vector<esglobal> rtoffset(threads);
//	#pragma omp parallel for
//	for (size_t t = 0; t < threads; t++) {
//		for (int r = rdistribution[t]; r < rdistribution[t + 1]; r++) {
//			roffset[r] += csize[r];
//			rtoffset[t] += roffset[r];
//		}
//	}
//
//	Esutils::sizesToOffsets(rtoffset);
//	#pragma omp parallel for
//	for (size_t t = 0; t < threads; t++) {
//		if (rdistribution[t] == rdistribution[t + 1]) {
//			continue;
//		}
//		esglobal tmp, prev = roffset[rdistribution[t]];
//		roffset[rdistribution[t]] = rtoffset[t];
//		for (int r = rdistribution[t] + 1; r < rdistribution[t + 1]; r++) {
//			tmp = roffset[r];
//			roffset[r] = roffset[r - 1] + prev;
//			prev = tmp;
//		}
//	}
//
//	MPI_Bcast(roffset.data(), environment->MPIsize * sizeof(esglobal), MPI_BYTE, environment->MPIsize - 1, MPI_COMM_WORLD);
//
//	#pragma omp parallel for
//	for (size_t t = 0; t < threads; t++) {
//		size_t offset;
//		for (int r = rdistribution[t]; r < rdistribution[t + 1]; r++) {
//			offset = 0;
//			for (size_t t = 0; t < threads; t++) {
//				offset += tsize[t][r];
//				tsize[t][r] = roffset[r] + coffset[r] + offset - tsize[t][r];
//			}
//		}
//	}
//
//	auto n2i = [&] (esglobal neighbor) {
//		return std::lower_bound(store->neighbors.begin(), store->neighbors.end(), neighbor) - store->neighbors.begin();
//	};
//
//	// threads x targets x (id, permutation)
//	std::vector<std::vector<std::vector<std::pair<esglobal, esglobal> > > > shalo(threads, std::vector<std::vector<std::pair<esglobal, esglobal> > >(store->neighbors.size()));
//	std::vector<std::vector<std::pair<esglobal, esglobal> > > rhalo(store->neighbors.size());
//
//	#pragma omp parallel for
//	for (size_t t = 0; t < threads; t++) {
//		std::vector<esglobal> toffset(tsize[t]);
//		for (size_t e = store->distribution[t]; e < store->distribution[t + 1]; e++) {
//			permutation[e] = toffset[partition[e]]++;
//		}
//
//		const crsiterator<esglobal> &dual = store->fullDual->iterator(t);
//		while (dual.next()) {
//			for (auto it = dual.begin(); it != dual.end(); ++it) {
//				if (*it < store->elementOffset || *it >= store->elementOffset + (esglobal)(store->distribution.back())) {
//					esglobal tindex = n2i(std::lower_bound(edistribution.begin(), edistribution.end(), *it + 1) - edistribution.begin() - 1);
//					shalo[t][tindex].push_back(std::make_pair(store->elementOffset + dual.index(), permutation[dual.index()]));
//				}
//			}
//		}
//	}
//
//	#pragma omp parallel for
//	for (size_t n = 0; n < store->neighbors.size(); n++) {
//		for (size_t t = 1; t < threads; t++) {
//			shalo[0][n].insert(shalo[0][n].end(), shalo[t][n].begin(), shalo[t][n].end());
//		}
//	}
//
//	if (!Communication::exchangeUnknownSize(shalo[0], rhalo, store->neighbors)) {
//		ESINFO(ERROR) << "ESPRESO internal error: exchanging dual halo elements.";
//	}
//
//	for (size_t n = 1; n < store->neighbors.size(); n++) {
//		rhalo[0].insert(rhalo[0].end(), rhalo[n].begin(), rhalo[n].end());
//	}
//
//	// threads x targets x data
//	std::vector<std::vector<std::vector<esglobal> > > tframes(threads, std::vector<std::vector<esglobal> >(targets.size()));
//	std::vector<std::vector<std::vector<esglobal> > > tneighs(threads, std::vector<std::vector<esglobal> >(targets.size()));
//	std::vector<std::vector<int> > tneighClusters(threads);
//
//	#pragma omp parallel for
//	for (size_t t = 0; t < threads; t++) {
//		const crsiterator<esglobal> &dual = store->fullDual->iterator(t);
//		while (dual.next()) {
//			size_t tindex = std::lower_bound(targets.begin(), targets.end(), partition[dual.index()]) - targets.begin();
//			tframes[t][tindex].push_back(dual.size());
//			for (auto it = dual.begin(); it != dual.end(); ++it) {
//				if (store->elementOffset <= *it && *it < store->elementOffset + (esglobal)(store->distribution.back())) {
//					tneighs[t][tindex].push_back(permutation[*it -store->elementOffset]);
//					tneighClusters[t].push_back(std::lower_bound(roffset.begin(), roffset.end(), tneighs[t][tindex].back() + 1) - roffset.begin() - 1);
//				} else {
//					std::pair<esglobal, esglobal> h(*it, 0);
//					auto halo = std::lower_bound(rhalo[0].begin(), rhalo[0].end(), h);
//					if (halo == rhalo[0].end()) {
//						ESINFO(ERROR) << "ESPRESO internal error: unknown halo element.";
//					}
//					tneighs[t][tindex].push_back(halo->second);
//				}
//			}
//		}
//		std::sort(tneighClusters[t].begin(), tneighClusters[t].end());
//		Esutils::removeDuplicity(tneighClusters[t]);
//	}
//
//	for (size_t t = 1; t < threads; t++) {
//		tneighClusters[0].insert(tneighClusters[0].end(), tneighClusters[t].begin(), tneighClusters[t].end());
//	}
//	std::sort(tneighClusters[0].begin(), tneighClusters[0].end());
//	Esutils::removeDuplicity(tneighClusters[0]);
//
//	#pragma omp parallel for
//	for (size_t n = 0; n < targets.size(); n++) {
//		for (size_t t = 1; t < threads; t++) {
//			tframes[0][n].insert(tframes[0][n].end(), tframes[t][n].begin(), tframes[t][n].end());
//			tneighs[0][n].insert(tneighs[0][n].end(), tneighs[t][n].begin(), tneighs[t][n].end());
//		}
//	}
//
//	std::vector<std::vector<esglobal> > rframes;
//	std::vector<std::vector<esglobal> > rneighs;
//
//	if (!(Communication::sendVariousTargets(tframes[0], rframes, targets) && Communication::sendVariousTargets(tneighs[0], rneighs, targets))) {
//		ESINFO(ERROR) << "ESPRESO internal error: distribute dual graph.";
//	}
//
//	for (size_t n = 1; n < rframes.size(); n++) {
//		rframes[0].insert(rframes[0].end(), rframes[n].begin(), rframes[n].end());
//		rneighs[0].insert(rneighs[0].end(), rneighs[n].begin(), rneighs[n].end());
//	}
//
//	if (!rframes.size()) {
//		// ParMETIS sometimes return empty cluster
//		rframes.resize(1);
//	}
//	std::vector<size_t> distribution = Esutils::getDistribution(threads, rframes[0].size());
//	store->distribution = distribution;
//
//	rframes[0].push_back(0);
//	for (size_t t = 0; t <= threads; t++) {
//		if (distribution[t] != 0 && distribution[t] == distribution[threads]) {
//			distribution[t]++;
//		}
//	}
//	parray<esglobal> frames(distribution.size(), distribution.data());
//
//	std::vector<esglobal> offsets(threads);
//	#pragma omp parallel for
//	for (size_t t = 0; t < threads; t++) {
//		esglobal offset = 0;
//		for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {
//			offset += rframes[0][e];
//		}
//		offsets[t] = offset;
//	}
//
//	Esutils::sizesToOffsets(offsets);
//
//	#pragma omp parallel for
//	for (size_t t = 0; t < threads; t++) {
//		esglobal offset = offsets[t];
//		for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {
//			frames[e] = offset;
//			offset += rframes[0][e];
//		}
//	}
//
//	for (size_t t = 0; t <= threads && frames.size(); t++) {
//		if (distribution[t] == frames.size()) {
//			distribution[t] = frames.back();
//		} else {
//			distribution[t] = frames[distribution[t]];
//		}
//	}
//
//	parray<esglobal> neighs(distribution.size(), distribution.data());
//
//	#pragma omp parallel for
//	for (size_t t = 0; t < threads; t++) {
//		if (distribution[t + 1] - distribution[t] > 0) {
//			memcpy(neighs.data() + distribution[t], rneighs[0].data() + distribution[t], sizeof(esglobal) * (distribution[t + 1] - distribution[t]));
//		}
//	}
//
//	if (mesh.__elements->fullDual != store->fullDual) {
//		delete store->fullDual;
//	}
//	store->fullDual = new pcrsdata<esglobal>(std::move(frames), std::move(neighs));
//	store->elementOffset = roffset[environment->MPIrank];
//
//	ESINFO(TVERBOSITY) << "Transformation::  distribution of Dual Graph finished.";
}

bool Transformation::_checkContinuity(NewMesh &mesh, ElementStore *store)
{
	ESINFO(TVERBOSITY) << "Transformation::  checking Continuity started.";


	ESINFO(TVERBOSITY) << "Transformation::  checking Continuity finished.";
	return true;
}

void Transformation::_tryrepartition(NewMesh &mesh, esglobal *permutation)
{
	ESINFO(TVERBOSITY) << "Transformation::  repartition of non-continuous partition started.";


	ESINFO(TVERBOSITY) << "Transformation::  repartition of non-continuous partition finished.";
}

void Transformation::_distributeNewMesh(NewMesh &mesh, ElementStore *store, esglobal *permutation)
{
//	ESINFO(TVERBOSITY) << "Transformation::  distribute mesh according new partition stated.";
//
//	std::vector<esglobal> edistribution(environment->MPIsize + 1);
//	MPI_Allgather(&store->elementOffset, sizeof(esglobal), MPI_BYTE, edistribution.data(), sizeof(esglobal), MPI_BYTE, MPI_COMM_WORLD);
//	edistribution.back() = store->elementOffset + store->distribution.back();
//	MPI_Bcast(&edistribution.back(), sizeof(esglobal), MPI_BYTE, environment->MPIsize - 1, MPI_COMM_WORLD);
//
//	if (store->intervals.size() < 1) {
//		store->intervals.push_back(std::vector<size_t>(edistribution.begin(), edistribution.end()));
//	} else {
//		store->intervals[0] = std::vector<size_t>(edistribution.begin(), edistribution.end());
//	}
//
//	size_t threads = environment->OMP_NUM_THREADS;
//
//	std::vector<std::vector<int> > ttargets(threads);
//	std::vector<int> targets;
//
//	#pragma omp parallel for
//	for (size_t t = 0; t < threads; t++) {
//		for (size_t e = mesh.__elements->distribution[t]; e < mesh.__elements->distribution[t + 1]; e++) {
//			if (permutation[e] < store->elementOffset || permutation[e] >= store->elementOffset + (esglobal)store->distribution.back()) {
//				ttargets[t].push_back(std::lower_bound(edistribution.begin(), edistribution.end(), permutation[e] + 1) - edistribution.begin() - 1);
//			}
//		}
//	}
//
//	for (size_t t = 0; t < threads; t++) {
//		targets.insert(targets.end(), ttargets[t].begin(), ttargets[t].end());
//	}
//	targets.push_back(environment->MPIrank);
//	std::sort(targets.begin(), targets.end());
//	Esutils::removeDuplicity(targets);
//
//	// threads x neighbors x data (number, node1, node2, node_number)
//	std::vector<std::vector<std::vector<esglobal> > > sBuffer(threads, std::vector<std::vector<esglobal> >(targets.size()));
//	std::vector<std::vector<esglobal> > rBuffer(targets.size());
//
//	#pragma omp parallel for
//	for (size_t t = 0; t < threads; t++) {
//		const crsiterator<eslocal> &indices = mesh.__elements->nodesIndices->iterator(t);
//
//		while(indices.next()) {
//			int target = std::lower_bound(edistribution.begin(), edistribution.end(), permutation[indices.index()] + 1) - edistribution.begin() - 1;
//			target = std::lower_bound(targets.begin(), targets.end(), target) - targets.begin();
//			sBuffer[t][target].push_back(-(*mesh.__elements->elements)[indices.index()]->vtkCode());
//			for (auto n = indices.begin(); n != indices.end(); ++n) {
//				sBuffer[t][target].push_back(mesh.coordinates().globalIndex(*n));
//			}
//		}
//	}
//
//	#pragma omp parallel for
//	for (size_t n = 0; n < targets.size(); n++) {
//		for (size_t t = 1; t < threads; t++) {
//			sBuffer[0][n].insert(sBuffer[0][n].end(), sBuffer[t][n].begin(), sBuffer[t][n].end());
//		}
//	}
//
//	if (!Communication::sendVariousTargets(sBuffer[0], rBuffer, targets)) {
//		ESINFO(ERROR) << "ESPRESO internal error: distribute mesh according new partition.";
//	}
//
//	for (size_t n = 1; n < rBuffer.size(); n++) {
//		rBuffer[0].insert(rBuffer[0].end(), rBuffer[n].begin(), rBuffer[n].end());
//	}
//
//	std::vector<size_t> rdistribution = Esutils::getDistribution(threads, rBuffer[0].size());
//	std::vector<std::vector<Element*> > elements(threads);
//
//	eslocal params[6] = { 0, 0, 0, 0, 0 };
//
//	#pragma omp parallel for
//	for (size_t t = 0; t < threads; t++) {
//		size_t p = rdistribution[t];
//		while (p < rdistribution[t + 1] && rBuffer[0][p] >= 0) {
//			p++;
//		}
//		while (p < rdistribution[t + 1]) {
//			switch (-rBuffer[0][p++]) {
//			case Line2VTKCode:
//				elements[t].push_back(new Line2(&rBuffer[0][p]));
//				p += Line2NodesCount;
//				break;
//			case Line3VTKCode:
//				elements[t].push_back(new Line3(&rBuffer[0][p]));
//				p += Line3NodesCount;
//				break;
//
//			case Triangle3VTKCode:
//				elements[t].push_back(new Triangle3(&rBuffer[0][p], params));
//				p += Triangle3NodesCount;
//				break;
//			case Triangle6VTKCode:
//				elements[t].push_back(new Triangle6(&rBuffer[0][p], params));
//				p += Triangle6NodesCount;
//				break;
//			case Square4VTKCode:
//				elements[t].push_back(new Square4(&rBuffer[0][p], params));
//				p += Square4NodesCount;
//				break;
//			case Square8VTKCode:
//				elements[t].push_back(new Square8(&rBuffer[0][p], params));
//				p += Square8NodesCount;
//				break;
//
//			case Hexahedron8VTKCode:
//				elements[t].push_back(new Hexahedron8(&rBuffer[0][p], Hexahedron8NodesCount, params));
//				p += Hexahedron8NodesCount;
//				break;
//			case Hexahedron20VTKCode:
//				elements[t].push_back(new Hexahedron20(&rBuffer[0][p], Hexahedron20NodesCount, params));
//				p += Hexahedron20NodesCount;
//				break;
//			case Tetrahedron4VTKCode:
//				elements[t].push_back(new Tetrahedron4(&rBuffer[0][p], Tetrahedron4NodesCount, params));
//				p += Tetrahedron4NodesCount;
//				break;
//			case Tetrahedron10VTKCode:
//				elements[t].push_back(new Tetrahedron10(&rBuffer[0][p], Tetrahedron10NodesCount, params));
//				p += Tetrahedron10NodesCount;
//				break;
//			case Prisma6VTKCode:
//				elements[t].push_back(new Prisma6(&rBuffer[0][p], Prisma6NodesCount, params));
//				p += Prisma6NodesCount;
//				break;
//			case Prisma15VTKCode:
//				elements[t].push_back(new Prisma15(&rBuffer[0][p], Prisma15NodesCount, params));
//				p += Prisma15NodesCount;
//				break;
//			case Pyramid5VTKCode:
//				elements[t].push_back(new Pyramid5(&rBuffer[0][p], Pyramid5NodesCount, params));
//				p += Pyramid5NodesCount;
//				break;
//			case Pyramid13VTKCode:
//				elements[t].push_back(new Pyramid13(&rBuffer[0][p], Pyramid13NodesCount, params));
//				p += Pyramid13NodesCount;
//				break;
//
//			default:
//				ESINFO(GLOBAL_ERROR) << "Unknown exchanged element.";
//			}
//		}
//	}
//
//	for (size_t t = 1; t < threads; t++) {
//		elements[0].insert(elements[0].end(), elements[t].begin(), elements[t].end());
//	}
//
//	store->elements = new parray<Element*>(store->distribution.size(), store->distribution.data());
//	#pragma omp parallel for
//	for (size_t t = 0; t < threads; t++) {
//		if (store->distribution[t + 1] - store->distribution[t] > 0) {
//			memcpy(store->elements->data() + store->distribution[t], elements[0].data() + store->distribution[t], sizeof(Element*) * (store->distribution[t + 1] - store->distribution[t]));
//		}
//	}
//
//	/// exchange permutation of neighbors elements
//
//	// threads x neighbors x
//	std::vector<std::vector<std::vector<std::pair<esglobal, esglobal> > > > sPermutation(threads, std::vector<std::vector<std::pair<esglobal, esglobal> > >(mesh.__elements->neighbors.size()));
//	std::vector<std::vector<std::pair<esglobal, esglobal> > > rPermutation(mesh.__elements->neighbors.size());
//
//
//	#pragma omp parallel for
//	for (size_t t = 0; t < threads; t++) {
//		const crsiterator<eslocal> &enodes = mesh.__elements->nodesIndices->iterator(t);
//		const crsiterator<eslocal> nranks = mesh.__nodes->ranks->iterator();
//
//		size_t first = mesh.__elements->elementOffset;
//		size_t last = mesh.__elements->elementOffset + mesh.__elements->elements->size();
//		const std::vector<int> &neighbors = mesh.__elements->neighbors;
//		std::vector<int> eranks;
//
//		while(enodes.next()) {
//			eranks.clear();
//			for (auto n = enodes.begin(); n != enodes.end(); ++n) {
//				nranks.seek(*n);
//				eranks.insert(eranks.end(), nranks.begin(), nranks.end());
//			}
//			std::sort(eranks.begin(), eranks.end());
//			Esutils::removeDuplicity(eranks);
//
//			for (auto r = eranks.begin(); r != eranks.end(); r++) {
//				if (*r != environment->MPIrank) {
//					size_t nindex = std::lower_bound(neighbors.begin(), neighbors.end(), *r) - neighbors.begin();
//					sPermutation[t][nindex].push_back(std::make_pair(enodes.index() + first, permutation[enodes.index()]));
//				}
//			}
//		}
//		for (size_t n = 0; n < neighbors.size(); n++) {
//			std::sort(sPermutation[t][n].begin(), sPermutation[t][n].end());
//			Esutils::removeDuplicity(sPermutation[t][n]);
//		}
//	}
//
//	#pragma omp parallel for
//	for (size_t n = 0; n < mesh.__elements->neighbors.size(); n++) {
//		for (size_t t = 1; t < threads; t++) {
//			sPermutation[0][n].insert(sPermutation[0][n].end(), sPermutation[t][n].begin(), sPermutation[t][n].end());
//		}
//		std::sort(sPermutation[0][n].begin(), sPermutation[0][n].end());
//		Esutils::removeDuplicity(sPermutation[0][n]);
//	}
//
//	if (!Communication::exchangeUnknownSize(sPermutation[0], rPermutation, mesh.__elements->neighbors)) {
//		ESINFO(ERROR) << "ESPRESO internal error: send permutation of neighbors elements.";
//	}
//
//	std::vector<esglobal> olddistribution(environment->MPIsize + 1);
//	MPI_Allgather(&mesh.__elements->elementOffset, sizeof(esglobal), MPI_BYTE, olddistribution.data(), sizeof(esglobal), MPI_BYTE, MPI_COMM_WORLD);
//	olddistribution.back() = mesh.__elements->elementOffset + mesh.__elements->distribution.back();
//	MPI_Bcast(&olddistribution.back(), sizeof(esglobal), MPI_BYTE, environment->MPIsize - 1, MPI_COMM_WORLD);
//
//	// update node parents vector
//
//	std::vector<std::vector<int> > neighbours(threads);
//	std::vector<std::vector<eslocal> > nranksCRS(threads);
//	std::vector<std::vector<eslocal> > nranksData(threads);
//	std::vector<std::vector<eslocal> > nParentsCRS(threads);
//	std::vector<std::vector<eslocal> > nParentsData(threads);
//
//	#pragma omp parallel for
//	for (size_t t = 0; t < threads; t++) {
//		const crsiterator<eslocal> &parents = mesh.__nodes->parentIndices->iterator(t);
//		const crsiterator<eslocal> oldranks = mesh.__nodes->ranks->iterator(t);
//		size_t first = mesh.__elements->elementOffset;
//		size_t last = mesh.__elements->elementOffset + mesh.__elements->elements->size();
//		const std::vector<int> &neighbors = mesh.__elements->neighbors;
//		std::pair<esglobal, esglobal> pair(0, 0);
//		std::vector<esglobal> nParents;
//		std::vector<eslocal> nranks;
//		esglobal eID;
//
//		while(parents.next() && oldranks.next()) {
//			nParents.clear();
//			nranks.clear();
//			for (auto it = parents.begin(); it != parents.end(); ++it) {
//				if (first <= *it && *it < last) {
//					eID = permutation[*it - first];
//				} else {
//					size_t cid = std::lower_bound(olddistribution.begin(), olddistribution.end(), *it + 1) - olddistribution.begin() - 1;
//					size_t nindex = std::lower_bound(neighbors.begin(), neighbors.end(), cid) - neighbors.begin();
//					pair.first = *it;
//					eID = std::lower_bound(rPermutation[nindex].begin(), rPermutation[nindex].end(), pair)->second;
//				}
//				nParents.push_back(eID);
//				size_t cid = std::lower_bound(edistribution.begin(), edistribution.end(), eID + 1) - edistribution.begin() - 1;
//				nranks.push_back(cid);
//			}
//
//			std::sort(nParents.begin(), nParents.end());
//			std::sort(nranks.begin(), nranks.end());
//			Esutils::removeDuplicity(nranks);
//
//			neighbours[t].insert(neighbours[t].end(), nranks.begin(), nranks.end());
//			std::sort(neighbours[t].begin(), neighbours[t].end());
//			Esutils::removeDuplicity(neighbours[t]);
//
//			if (oldranks.front() == environment->MPIrank) {
//				nranksCRS[t].push_back(nranks.size());
//				nranksData[t].insert(nranksData[t].end(), nranks.begin(), nranks.end());
//
//				nParentsCRS[t].push_back(nParents.size());
//				nParentsData[t].insert(nParentsData[t].end(), nParents.begin(), nParents.end());
//			}
//		}
//	}
//
//	for (size_t t = 1; t < threads; t++) {
//		neighbours[0].insert(neighbours[0].end(), neighbours[t].begin(), neighbours[t].end());
//	}
//	std::sort(neighbours[0].begin(), neighbours[0].end());
//	Esutils::removeDuplicity(neighbours[0]);
//
//	std::vector<std::vector<std::vector<esglobal> > > sNodeIDs(threads, std::vector<std::vector<esglobal> >(neighbours[0].size()));
//	std::vector<std::vector<std::vector<esglobal> > > sNodeParentsCRS(threads, std::vector<std::vector<esglobal> >(neighbours[0].size()));
//	std::vector<std::vector<std::vector<esglobal> > > sNodeParentsData(threads, std::vector<std::vector<esglobal> >(neighbours[0].size()));
//	std::vector<std::vector<std::vector<double> > > sNodeCoordinates(threads, std::vector<std::vector<double> >(neighbours[0].size()));
//
//	std::vector<std::vector<esglobal> > rNodeIDs;
//	std::vector<std::vector<esglobal> > rNodeParentsCRS;
//	std::vector<std::vector<esglobal> > rNodeParentsData;
//	std::vector<std::vector<double> > rNodeCoordinates;
//
//	#pragma omp parallel for
//	for (size_t t = 0; t < threads; t++) {
//		const crsiterator<eslocal> oldranks = mesh.__nodes->ranks->iterator(t);
//		size_t pRank = 0, pRankData = 0;
//		size_t pParent = 0, pParentData = 0;
//
//		while(oldranks.next()) {
//			if (oldranks.front() == environment->MPIrank) {
//				for (size_t r = 0; r < nranksCRS[t][pRank]; r++) {
//					size_t tindex = std::lower_bound(neighbours[0].begin(), neighbours[0].end(), nranksData[t][pRankData + r]) - neighbours[0].begin();
//					sNodeIDs[t][tindex].push_back(mesh.coordinates().globalIndex(oldranks.index()));
//					sNodeParentsCRS[t][tindex].push_back(nParentsCRS[t][pParent]);
//					sNodeParentsData[t][tindex].insert(sNodeParentsData[t][tindex].end(), &nParentsData[t][pParentData], &nParentsData[t][pParentData + nParentsCRS[t][pParent]]);
//					sNodeCoordinates[t][tindex].push_back(mesh.coordinates()[oldranks.index()].x);
//					sNodeCoordinates[t][tindex].push_back(mesh.coordinates()[oldranks.index()].y);
//					if (mesh.__elements->dimension == 3) {
//						sNodeCoordinates[t][tindex].push_back(mesh.coordinates()[oldranks.index()].z);
//					}
//				}
//				pRankData += nranksCRS[t][pRank++];
//				pParentData += nParentsCRS[t][pParent++];
//			}
//		}
//	}
//
//	#pragma omp parallel for
//	for (size_t n = 0; n < neighbours[0].size(); n++) {
//		for (size_t t = 1; t < threads; t++) {
//			sNodeIDs[0][n].insert(sNodeIDs[0][n].end(), sNodeIDs[t][n].begin(), sNodeIDs[t][n].end());
//			sNodeParentsCRS[0][n].insert(sNodeParentsCRS[0][n].end(), sNodeParentsCRS[t][n].begin(), sNodeParentsCRS[t][n].end());
//			sNodeParentsData[0][n].insert(sNodeParentsData[0][n].end(), sNodeParentsData[t][n].begin(), sNodeParentsData[t][n].end());
//			sNodeCoordinates[0][n].insert(sNodeCoordinates[0][n].end(), sNodeCoordinates[t][n].begin(), sNodeCoordinates[t][n].end());
//		}
//	}
//
//	if (!(
//			Communication::sendVariousTargets(sNodeIDs[0], rNodeIDs, neighbours[0]) &&
//			Communication::sendVariousTargets(sNodeParentsCRS[0], rNodeParentsCRS, neighbours[0]) &&
//			Communication::sendVariousTargets(sNodeParentsData[0], rNodeParentsData, neighbours[0]) &&
//			Communication::sendVariousTargets(sNodeCoordinates[0], rNodeCoordinates, neighbours[0])
//			)) {
//
//		ESINFO(ERROR) << "ESPRESO internal error: distribute nodes according new partition.";
//	}
//
//
//	std::vector<size_t> rOffsets(rNodeIDs.size()), rDataOffsets(rNodeIDs.size());
//	for (size_t n = 0; n < rNodeIDs.size(); n++) {
//		rOffsets[n] = rNodeIDs[n].size();
//		rDataOffsets[n] = rNodeParentsData[n].size();
//	}
//	size_t rSize = Esutils::sizesToOffsets(rOffsets);
//	Esutils::sizesToOffsets(rDataOffsets);
//	delete mesh.__nodes->ranks;
//	if (mesh.__nodes->globalIDs != NULL) {
//		delete mesh.__nodes->globalIDs;
//	}
//	if (mesh.__nodes->coordinates != NULL) {
//		delete mesh.__nodes->coordinates;
//	}
//	if (mesh.__nodes->parentIndices != NULL) {
//		delete mesh.__nodes->parentIndices;
//	}
//
//	mesh.__nodes->globalIDs = new parray<esglobal>(threads, rSize);
//	std::vector<size_t> cdistribution(mesh.__nodes->globalIDs->distribution(), mesh.__nodes->globalIDs->distribution() + threads + 1);
//	std::vector<size_t> pdistribution(mesh.__nodes->globalIDs->distribution(), mesh.__nodes->globalIDs->distribution() + threads + 1);
//	for (size_t t = 0; t < cdistribution.size(); t++) {
//		cdistribution[t] *= mesh.__nodes->dimension;
//		if (pdistribution[t] == pdistribution.back()) {
//			pdistribution[t]++;
//		}
//	}
//
//	size_t dimension = mesh.__nodes->dimension;
//	mesh.__nodes->coordinates = new parray<double>(cdistribution.size(), cdistribution.data());
//	parray<esglobal> eindices(pdistribution.size(), pdistribution.data());
//
//	#pragma omp parallel for
//	for (size_t t = 0; t < threads; t++) {
//		if (t + 1 == threads) {
//			eindices.back() = 0;
//		}
//		if (mesh.__nodes->globalIDs->distribution()[t + 1] == 0) {
//			continue;
//		}
//		size_t begin = mesh.__nodes->globalIDs->distribution()[t];
//		size_t end   = mesh.__nodes->globalIDs->distribution()[t + 1];
//
//		size_t beginsource = std::lower_bound(rOffsets.begin(), rOffsets.end(), begin + 1) - rOffsets.begin() - 1;
//		size_t endsource   = std::upper_bound(rOffsets.begin(), rOffsets.end(), end) - rOffsets.begin() - 1;
//
//		size_t tstart = begin - rOffsets[beginsource];
//		size_t tend   = end   - rOffsets[endsource];
//		size_t tsize  = beginsource == endsource && rNodeIDs[beginsource].size() > tend ? tend - tstart : rNodeIDs[beginsource].size() - tstart;
//		memcpy(mesh.__nodes->globalIDs->data() + begin, rNodeIDs[beginsource].data() + tstart, sizeof(esglobal) * tsize);
//		memcpy(mesh.__nodes->coordinates->data() + dimension * begin, rNodeCoordinates[beginsource].data() + dimension * tstart, dimension * sizeof(double) * tsize);
//		memcpy(eindices.data() + begin, rNodeParentsCRS[beginsource].data() + tstart, sizeof(esglobal) * tsize);
//		for (size_t i = beginsource + 1; i < endsource; i++) {
//			memcpy(mesh.__nodes->globalIDs->data() + rOffsets[i], rNodeIDs[i].data(), sizeof(esglobal) * rNodeIDs[i].size());
//			memcpy(mesh.__nodes->coordinates->data() + dimension * rOffsets[i], rNodeCoordinates[i].data(), dimension * sizeof(double) * rNodeIDs[i].size());
//			memcpy(eindices.data() + rOffsets[i], rNodeParentsCRS[i].data(), sizeof(esglobal) * rNodeParentsCRS[i].size());
//		}
//		if (beginsource != endsource) {
//			memcpy(mesh.__nodes->globalIDs->data() + rOffsets[endsource], rNodeIDs[endsource].data(), sizeof(esglobal) * (end - rOffsets[endsource]));
//			memcpy(mesh.__nodes->coordinates->data() + dimension * rOffsets[endsource], rNodeCoordinates[endsource].data(), dimension * sizeof(double) * (end - rOffsets[endsource]));
//			memcpy(eindices.data() + rOffsets[endsource], rNodeParentsCRS[endsource].data(), sizeof(esglobal) * (end - rOffsets[endsource]));
//		}
//	}
//
//
//	std::vector<size_t> eindicesOffsets(threads);
//	#pragma omp parallel for
//	for (size_t t = 0; t < threads; t++) {
//		size_t esize = 0;
//		for (size_t e = eindices.distribution()[t]; e < eindices.distribution()[t + 1]; e++) {
//			esize += eindices[e];
//		}
//		eindicesOffsets[t] = esize;
//	}
//	Esutils::sizesToOffsets(eindicesOffsets);
//	#pragma omp parallel for
//	for (size_t t = 0; t < threads; t++) {
//		size_t eoffset = eindicesOffsets[t];
//		for (size_t e = eindices.distribution()[t]; e < eindices.distribution()[t + 1]; e++) {
//			eoffset += eindices[e];
//			eindices[e] = eoffset - eindices[e];
//		}
//	}
//
//	std::vector<size_t> datadistribution(pdistribution);
//	for (size_t t = 0; t < datadistribution.size(); t++) {
//		if (datadistribution[t] == datadistribution.back()) {
//			datadistribution[t] = eindices.back();
//		} else {
//			datadistribution[t] = eindices[datadistribution[t]];
//		}
//
//	}
//
//	parray<esglobal> eindicesdata(datadistribution.size(), datadistribution.data());
//	#pragma omp parallel for
//	for (size_t t = 0; t < threads; t++) {
//		if (datadistribution[t + 1] == 0) {
//			continue;
//		}
//		size_t begin = datadistribution[t];
//		size_t end   = datadistribution[t + 1];
//
//		size_t beginsource = std::lower_bound(rDataOffsets.begin(), rDataOffsets.end(), begin + 1) - rDataOffsets.begin() - 1;
//		size_t endsource   = std::upper_bound(rDataOffsets.begin(), rDataOffsets.end(), end) - rDataOffsets.begin() - 1;
//
//		size_t tstart = begin - rDataOffsets[beginsource];
//		size_t tend   = end   - rDataOffsets[endsource];
//		size_t tsize  = beginsource == endsource && rNodeParentsData[beginsource].size() > tend ? tend - tstart : rNodeParentsData[beginsource].size() - tstart;
//		memcpy(eindicesdata.data() + begin, rNodeParentsData[beginsource].data() + tstart, sizeof(esglobal) * tsize);
//		for (size_t i = beginsource + 1; i < endsource; i++) {
//			memcpy(eindicesdata.data() + rDataOffsets[i], rNodeParentsData[i].data(), sizeof(esglobal) * rNodeParentsData[i].size());
//		}
//		if (beginsource != endsource) {
//			memcpy(eindicesdata.data() + rDataOffsets[endsource], rNodeParentsData[endsource].data(), sizeof(esglobal) * (end - rDataOffsets[endsource]));
//		}
//	}
//
//	mesh.__nodes->parentIndices = new pcrsdata<esglobal>(std::move(eindices), std::move(eindicesdata));
//
//	eslocal *npermutation = new eslocal[mesh.__nodes->globalIDs->size()];
//	std::iota(npermutation, npermutation + mesh.__nodes->globalIDs->size(), 0);
//	std::sort(npermutation, npermutation + mesh.__nodes->globalIDs->size(), [&] (eslocal i, eslocal j) { return (*mesh.__nodes->globalIDs)[i] < (*mesh.__nodes->globalIDs)[j]; });
//
//	mesh.__nodes->globalIDs->permute(npermutation);
//	mesh.__nodes->parentIndices->permute(npermutation);
//
//	parray<double> pcenters(mesh.__nodes->coordinates->threads(), mesh.__nodes->coordinates->distribution());
//	#pragma omp parallel for
//	for (size_t t = 0; t < threads; t++) {
//		for (size_t i = mesh.__nodes->coordinates->distribution()[t] / dimension; i < mesh.__nodes->coordinates->distribution()[t + 1] / dimension; i++) {
//			pcenters[dimension * i] = (*mesh.__nodes->coordinates)[dimension * npermutation[i]];
//			pcenters[dimension * i + 1] = (*mesh.__nodes->coordinates)[dimension * npermutation[i] + 1];
//			if (dimension == 3) {
//				pcenters[dimension * i + 2] = (*mesh.__nodes->coordinates)[dimension * npermutation[i] + 2];
//			}
//		}
//	}
//	delete[] npermutation;
//
//	mesh.__nodes->coordinates->swap(pcenters);
//
//	std::vector<std::vector<eslocal> > indicesSize(threads);
//	std::vector<std::vector<eslocal> > indices(threads);
//	#pragma omp parallel for
//	for (size_t t = 0; t < threads; t++) {
//		for (size_t e = store->elements->distribution()[t]; e < store->elements->distribution()[t + 1]; e++) {
//			indicesSize[t].push_back((*store->elements)[e]->nodes());
//			for (size_t n = 0; n < (*store->elements)[e]->nodes(); n++) {
//				eslocal index = std::lower_bound(mesh.__nodes->globalIDs->begin(), mesh.__nodes->globalIDs->end(), (*store->elements)[e]->node(n)) - mesh.__nodes->globalIDs->begin();
//				(*store->elements)[e]->node(n) = index;
//				indices[t].push_back(index);
//			}
//		}
//	}
//
//	store->nodesIndices = new pcrsdata<eslocal>(indicesSize, indices);
//
//	std::vector<std::vector<eslocal> > rankSize(threads);
//	std::vector<std::vector<eslocal> > ranks(threads);
//	std::vector<std::vector<eslocal> > tranks(threads);
//	#pragma omp parallel for
//	for (size_t t = 0; t < threads; t++) {
//		const crsiterator<esglobal> &eit = mesh.__nodes->parentIndices->iterator(t);
//		std::vector<eslocal> eranks;
//		while (eit.next()) {
//			eranks.clear();
//			for (auto e = eit.begin(); e != eit.end(); ++e) {
//				eranks.push_back(std::lower_bound(edistribution.begin(), edistribution.end(), *e + 1) - edistribution.begin() - 1);
//			}
//			std::sort(eranks.begin(), eranks.end());
//			Esutils::removeDuplicity(eranks);
//			rankSize[t].push_back(eranks.size());
//			ranks[t].insert(ranks[t].end(), eranks.begin(), eranks.end());
//
//			tranks[t].insert(tranks[t].end(), eranks.begin(), eranks.end());
//			std::sort(tranks[t].begin(), tranks[t].end());
//			Esutils::removeDuplicity(tranks[t]);
//		}
//	}
//
//	mesh.__nodes->ranks = new pcrsdata<eslocal>(rankSize, ranks);
//	for (size_t t = 1; t < threads; t++) {
//		tranks[0].insert(tranks[0].end(), tranks[t].begin(), tranks[t].end());
//	}
//	std::sort(tranks[0].begin(), tranks[0].end());
//	Esutils::removeDuplicity(tranks[0]);
//	store->neighbors.clear();
//	for (size_t i = 0; i < tranks[0].size(); i++) {
//		if (tranks[0][i] != environment->MPIrank) {
//			store->neighbors.push_back(tranks[0][i]);
//		}
//	}
//
//	mesh.__nodes->distribution = std::vector<size_t>(mesh.__nodes->globalIDs->distribution(), mesh.__nodes->globalIDs->distribution() + mesh.__nodes->globalIDs->threads() + 1);
//
//	ESINFO(TVERBOSITY) << "Transformation::  distribute mesh according new partition finished.";
}

void Transformation::_permuteNewMesh(NewMesh &mesh, esglobal *permutation)
{
//	ESINFO(TVERBOSITY) << "Transformation::  permutation of mesh elements started.";
//
//	size_t threads = environment->OMP_NUM_THREADS;
//
//	esglobal *backpermutation = new esglobal[mesh.__elements->elements->size()];
//	for (size_t i = 0; i < mesh.__elements->elements->size(); i++) {
//		backpermutation[permutation[i]] = i;
//	}
//
//	if (mesh.__elements->fullDual != NULL || mesh.__nodes->parentIndices != NULL) {
//
//		std::vector<esglobal> edistribution(environment->MPIsize + 1);
//		MPI_Allgather(&mesh.__elements->elementOffset, sizeof(esglobal), MPI_BYTE, edistribution.data(), sizeof(esglobal), MPI_BYTE, MPI_COMM_WORLD);
//		edistribution.back() = mesh.__elements->elementOffset + mesh.__elements->elements->size();
//		MPI_Bcast(&edistribution.back(), sizeof(esglobal), MPI_BYTE, environment->MPIsize - 1, MPI_COMM_WORLD);
//
//		std::vector<std::vector<std::vector<std::pair<esglobal, esglobal> > > > sBuffer(threads, std::vector<std::vector<std::pair<esglobal, esglobal> > >(mesh.neighbours().size()));
//		std::vector<std::vector<std::pair<esglobal, esglobal> > > rBuffer(mesh.neighbours().size());
//
//		esglobal eoffset = mesh.__elements->elementOffset;
//
//		if (mesh.__nodes->parentIndices != NULL) {
//			#pragma omp parallel for
//			for (size_t t = 0; t < threads; t++) {
//				const crsiterator<esglobal> &element = mesh.__elements->nodesIndices->iterator(t);
//				const crsiterator<eslocal> nranks = mesh.__nodes->ranks->iterator();
//				size_t neighbor;
//				while (element.next()) {
//					for (auto n = element.begin(); n != element.end(); ++n) {
//						nranks.seek(*n);
//						if (nranks.size() > 1) {
//							for (auto c = nranks.begin(); c != nranks.end(); ++c) {
//								if (*c != environment->MPIrank) {
//									neighbor = std::lower_bound(mesh.neighbours().begin(), mesh.neighbours().end(), *c) - mesh.neighbours().begin();
//									sBuffer[t][neighbor].push_back(std::make_pair(eoffset + element.index(), backpermutation[element.index()] + eoffset));
//								}
//							}
//						}
//					}
//				}
//			}
//		} else {
//			#pragma omp parallel for
//			for (size_t t = 0; t < threads; t++) {
//				const crsiterator<esglobal> &element = mesh.__elements->fullDual->iterator(t);
//				size_t neighbor;
//				while (element.next()) {
//					for (auto n = element.begin(); n != element.end(); ++n) {
//						if (*n - eoffset >= mesh.__elements->elements->size()) {
//							neighbor = std::lower_bound(edistribution.begin(), edistribution.end(), *n + 1) - edistribution.begin() - 1;
//							neighbor = std::lower_bound(mesh.neighbours().begin(), mesh.neighbours().end(), neighbor) - mesh.neighbours().begin();
//							sBuffer[t][neighbor].push_back(std::make_pair(eoffset + element.index(), backpermutation[element.index()] + eoffset));
//						}
//					}
//				}
//			}
//		}
//
//		#pragma omp parallel for
//		for (size_t t = 0; t < threads; t++) {
//			for (size_t n = 0; n < mesh.neighbours().size(); n++) {
//				std::sort(sBuffer[t][n].begin(), sBuffer[t][n].end());
//				Esutils::removeDuplicity(sBuffer[t][n]);
//			}
//		}
//
//		#pragma omp parallel for
//		for (size_t n = 0; n < mesh.neighbours().size(); n++) {
//			for (size_t t = 1; t < threads; t++) {
//				sBuffer[0][n].insert(sBuffer[0][n].end(), sBuffer[t][n].begin(), sBuffer[t][n].end());
//			}
//			std::sort(sBuffer[0][n].begin(), sBuffer[0][n].end());
//			Esutils::removeDuplicity(sBuffer[0][n]);
//		}
//
//		if (!Communication::exchangeUnknownSize(sBuffer[0], rBuffer, mesh.neighbours())) {
//			ESINFO(ERROR) << "ESPRESO internal error while exchanging mesh permutation.";
//		}
//
//		std::vector<pcrsdata<esglobal>* > items;
//		if (mesh.__elements->fullDual != NULL) {
//			items.push_back(mesh.__elements->fullDual);
//			mesh.__elements->fullDual->permute(permutation);
//		}
//		if (mesh.__nodes->parentIndices != NULL) {
//			items.push_back(mesh.__nodes->parentIndices);
//		}
//
//		for (size_t i = 0; i < items.size(); i++) {
//			#pragma omp parallel for
//			for (size_t t = 0; t < threads; t++) {
//				crsiterator<esglobal> &item = items[i]->iterator(t);
//				size_t neighbor;
//				std::pair<esglobal, esglobal> tmpPair(0, -1);
//				while (item.next()) {
//					for (auto n = item.begin(); n != item.end(); ++n) {
//						if (*n - eoffset >= mesh.__elements->elements->size()) {
//							neighbor = std::lower_bound(edistribution.begin(), edistribution.end(), *n + 1) - edistribution.begin() - 1;
//							neighbor = std::lower_bound(mesh.neighbours().begin(), mesh.neighbours().end(), neighbor) - mesh.neighbours().begin();
//							tmpPair.first = *n;
//							*n = std::lower_bound(rBuffer[neighbor].begin(), rBuffer[neighbor].end(), tmpPair)->second;
//						} else {
//							*n = backpermutation[*n - eoffset] + eoffset;
//						}
//					}
//				}
//			}
//		}
//
//		#pragma omp parallel for
//		for (size_t t = 0; t < threads; t++) {
//			crsiterator<esglobal> &rdual = mesh.__elements->restrictedDual->iterator(t);
//			while (rdual.next()) {
//				for (auto n = rdual.begin(); n != rdual.end(); ++n) {
//					*n = backpermutation[*n];
//				}
//			}
//		}
//	}
//
//	// WARNING: elements are created in different memory location (NUMA). re-allocation??
//	mesh.__elements->elements->permute(permutation);
//	mesh.__elements->nodesIndices->permute(permutation);
//	mesh.__elements->restrictedDual->permute(permutation);
//
//	delete[] backpermutation;
//	ESINFO(TVERBOSITY) << "Transformation::  permutation of mesh elements finished.";
}


