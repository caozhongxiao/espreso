
#include "transformations.h"

#include "../newmesh.h"

#include "../../basis/logging/logging.h"

using namespace espreso;

void Transformation::linkElementsFromTo(NewMesh &mesh, TFlags::ELEVEL from, TFlags::ELEVEL to)
{
	auto elevel = [] (TFlags::ELEVEL level) {
		switch (level) {
		case TFlags::ELEVEL::ELEMENT:
			return "ELEMENTS";
		case TFlags::ELEVEL::FACE:
			return "FACES";
		case TFlags::ELEVEL::EDGE:
			return "EDGES";
		case TFlags::ELEVEL::NODE:
			return "NODES";
		default:
			return "??";
		}
	};
	ESINFO(TVERBOSITY) << "Transformation::add link from " << elevel(from) << " to " << elevel(to) << " started\n";
//	ElementStore *pstore;
//	ElementStore *chstore;
//	switch (parents) {
//	case TFlags::ELEVEL::ELEMENT:
//		pstore = mesh.__elements;
//		break;
//	case TFlags::ELEVEL::FACE:
//		pstore = mesh.__faces;
//		break;
//	case TFlags::ELEVEL::EDGE:
//		pstore = mesh.__edges;
//		break;
//	default:
//		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: invalid parent elements.";
//	}
//
//	switch (children) {
//	case TFlags::ELEVEL::FACE:
//	case TFlags::ELEVEL::EDGE:
//		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: implement fillParentElements.";
//		break;
//	case TFlags::ELEVEL::NODE:
//		chstore = mesh.__nodes;
//		break;
//	default:
//		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: invalid children elements.";
//	}
//
//
//	if (!pstore->elements->size()) {
//		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: fillParentElements assumes non-empty parent vector.";
//	}
//
//	const parray<Element*> &elements = *pstore->elements;
//	size_t threads = environment->OMP_NUM_THREADS;
//
//	auto n2i = [ & ] (size_t neighbour) {
//		return std::lower_bound(mesh.neighbours().begin(), mesh.neighbours().end(), neighbour) - mesh.neighbours().begin();
//	};
//
//	// thread x neighbor x (node, element)
//	std::vector<std::vector<std::vector<std::pair<esglobal, esglobal> > > > sBuffer(threads, std::vector<std::vector<std::pair<esglobal, esglobal> > >(mesh.neighbours().size()));
//	std::vector<std::vector<std::pair<esglobal, esglobal> > > rBuffer(mesh.neighbours().size());
//	std::vector<std::vector<std::pair<esglobal, esglobal> > > nepairs(threads);
//
//	#pragma omp parallel for
//	for (size_t t = 0; t < threads; t++) {
//		const crsiterator<eslocal>& enodes = pstore->nodesIndices->iterator(t);
//		const crsiterator<eslocal> nranks = chstore->ranks->iterator();
//		while (enodes.next()) {
//			for (auto n = enodes.begin(); n != enodes.end(); ++n) {
//				nepairs[t].push_back(std::make_pair(mesh.coordinates().globalIndex(*n), enodes.index() + pstore->elementOffset));
//
//				nranks.seek(*n);
//				for (auto c = nranks.begin(); c != nranks.end(); ++c) {
//					if (*c != environment->MPIrank) {
//						sBuffer[t][n2i(*c)].push_back(std::make_pair(mesh.coordinates().globalIndex(*n), enodes.index() + pstore->elementOffset));
//					}
//				}
//			}
//		}
//		std::sort(nepairs[t].begin(), nepairs[t].end());
//	}
//
//	#pragma omp parallel for
//	for (size_t n = 0; n < sBuffer[0].size(); n++) {
//		for (size_t t = 1; t < threads; t++) {
//			sBuffer[0][n].insert(sBuffer[0][n].end(), sBuffer[t][n].begin(), sBuffer[t][n].end());
//		}
//		std::sort(sBuffer[0][n].begin(), sBuffer[0][n].end());
//	}
//
//	if (!Communication::exchangeUnknownSize(sBuffer[0], rBuffer, mesh.neighbours())) {
//		ESINFO(ERROR) << "ESPRESO internal error while exchanging parent elements.";
//	}
//
//	std::vector<std::vector<std::pair<esglobal, esglobal> > > nparents(threads);
//
//	#pragma omp parallel for
//	for (size_t t = 0; t < threads; t++) {
//		if (chstore->distribution[t] == chstore->distribution[t + 1]) {
//			continue;
//		}
//		std::pair<esglobal, esglobal> first(mesh.coordinates().globalIndex(chstore->distribution[t]), 0);
//		std::pair<esglobal, esglobal>  last(mesh.coordinates().globalIndex(chstore->distribution[t + 1] - 1), 0);
//
//		auto compare = [] (const std::pair<esglobal, esglobal> &p, const std::pair<esglobal, esglobal> &node) { return p.first < node.first; };
//
//		for (size_t n = 0; n < rBuffer.size(); n++) {
//			auto begin = std::lower_bound(rBuffer[n].begin(), rBuffer[n].end(), first, compare);
//			auto end   = std::upper_bound(rBuffer[n].begin(), rBuffer[n].end(), last , compare);
//			if (begin < end) {
//				nparents[t].insert(nparents[t].end(), begin, end);
//			}
//		}
//		for (size_t p = 0; p < threads; p++) {
//			auto begin = std::lower_bound(nepairs[p].begin(), nepairs[p].end(), first, compare);
//			auto end   = std::upper_bound(nepairs[p].begin(), nepairs[p].end(), last , compare);
//			if (begin < end) {
//				nparents[t].insert(nparents[t].end(), begin, end);
//			}
//		}
//		std::sort(nparents[t].begin(), nparents[t].end());
//	}
//
//	std::vector<size_t> odistribution = { 0 };
//	std::vector<size_t> edistribution = { 0 };
//	for (size_t t = 0; t < threads; t++) {
//		odistribution.push_back(odistribution.back() + chstore->distribution[t + 1] - chstore->distribution[t]);
//		edistribution.push_back(edistribution.back() + nparents[t].size());
//	}
//	odistribution.back() += 1;
//
//	parray<esglobal> offsets(odistribution.size(), odistribution.data());
//	parray<esglobal> nelements(edistribution.size(), odistribution.data());
//
//	#pragma omp parallel for
//	for (size_t t = 0; t < threads; t++) {
//		size_t index = chstore->distribution[t];
//		size_t offset = edistribution[t];
//
//		offsets[index] = offset;
//		for (size_t n = 0; n < nparents[t].size(); n++, offset++) {
//			nelements[offset] = nparents[t][n].second;
//			if (n && nparents[t][n - 1].first != nparents[t][n].first) {
//				offsets[++index] = offset;
//			}
//		}
//	}
//	offsets[offsets.size() - 1] = edistribution.back();
//
//	chstore->parentIndices = new pcrsdata<esglobal>(std::move(offsets), std::move(nelements));

	ESINFO(TVERBOSITY) << "Transformation::add link from " << elevel(from) << " to " << elevel(to) << " finished\n";
}

void Transformation::fillHaloElement(NewMesh &mesh)
{
//	const parray<Element*> &elements = *mesh.__elements->elements;
//	size_t threads = environment->OMP_NUM_THREADS;
//
//	auto n2i = [ & ] (size_t neighbour) {
//		return std::lower_bound(mesh.neighbours().begin(), mesh.neighbours().end(), neighbour) - mesh.neighbours().begin();
//	};
//
//	// thread x neighbor x (node, element)
//	std::vector<std::vector<std::vector<HaloElement> > > sBuffer(threads, std::vector<std::vector<HaloElement> >(mesh.neighbours().size()));
//	std::vector<std::vector<HaloElement> > rBuffer(mesh.neighbours().size());
//
//	#pragma omp parallel for
//	for (size_t t = 0; t < threads; t++) {
//		const crsiterator<eslocal> &enodes = mesh.__elements->nodesIndices->iterator(t);
//		const crsiterator<eslocal> nranks = mesh.__nodes->ranks->iterator();
//		HaloElement h;
//
//		while (enodes.next()) {
//			h.id       = enodes.index() + mesh.__elements->elementOffset;
//			h.nCommon  = elements[enodes.index()]->nCommon();
//			h.body     = elements[enodes.index()]->param(Element::BODY);
//			h.material = elements[enodes.index()]->param(Element::MATERIAL);
//			h.type     = (int)elements[enodes.index()]->type();
//			for (auto n = enodes.begin(); n != enodes.end(); ++n) {
//				nranks.seek(*n);
//				if (nranks.size() > 1) {
//					for (auto c = nranks.begin(); c != nranks.end(); ++c) {
//						if (*c != environment->MPIrank) {
//							sBuffer[t][n2i(*c)].push_back(h);
//						}
//					}
//				}
//			}
//		}
//		for (size_t n = 0; n < mesh.neighbours().size(); n++) {
//			std::sort(sBuffer[t][n].begin(), sBuffer[t][n].end());
//			Esutils::removeDuplicity(sBuffer[t][n]);
//		}
//	}
//
//	#pragma omp parallel for
//	for (size_t n = 0; n < sBuffer[0].size(); n++) {
//		for (size_t t = 1; t < threads; t++) {
//			sBuffer[0][n].insert(sBuffer[0][n].end(), sBuffer[t][n].begin(), sBuffer[t][n].end());
//		}
//		std::sort(sBuffer[0][n].begin(), sBuffer[0][n].end());
//		Esutils::removeDuplicity(sBuffer[0][n]);
//	}
//
//	if (!Communication::exchangeUnknownSize(sBuffer[0], rBuffer, mesh.neighbours())) {
//		ESINFO(ERROR) << "ESPRESO internal error while exchanging parent elements.";
//	}
//
//	mesh.__elements->haloElements = new std::vector<HaloElement>();
//	for (size_t n = 0; n < rBuffer.size(); n++) {
//		mesh.__elements->haloElements->insert(mesh.__elements->haloElements->end(), rBuffer[n].begin(), rBuffer[n].end());
//		std::sort(mesh.__elements->haloElements->begin(), mesh.__elements->haloElements->end());
//		Esutils::removeDuplicity(*mesh.__elements->haloElements);
//	}
}

void Transformation::fillParentElements(NewMesh &mesh, TFlags::ELEVEL children, TFlags::ELEVEL parents)
{
//	ElementStore *pstore;
//	ElementStore *chstore;
//	switch (parents) {
//	case TFlags::ELEVEL::ELEMENT:
//		pstore = mesh.__elements;
//		break;
//	case TFlags::ELEVEL::FACE:
//		pstore = mesh.__faces;
//		break;
//	case TFlags::ELEVEL::EDGE:
//		pstore = mesh.__edges;
//		break;
//	default:
//		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: invalid parent elements.";
//	}
//
//	switch (children) {
//	case TFlags::ELEVEL::FACE:
//	case TFlags::ELEVEL::EDGE:
//		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: implement fillParentElements.";
//		break;
//	case TFlags::ELEVEL::NODE:
//		chstore = mesh.__nodes;
//		break;
//	default:
//		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: invalid children elements.";
//	}
//
//
//	if (!pstore->elements->size()) {
//		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: fillParentElements assumes non-empty parent vector.";
//	}
//
//	const parray<Element*> &elements = *pstore->elements;
//	size_t threads = environment->OMP_NUM_THREADS;
//
//	auto n2i = [ & ] (size_t neighbour) {
//		return std::lower_bound(mesh.neighbours().begin(), mesh.neighbours().end(), neighbour) - mesh.neighbours().begin();
//	};
//
//	// thread x neighbor x (node, element)
//	std::vector<std::vector<std::vector<std::pair<esglobal, esglobal> > > > sBuffer(threads, std::vector<std::vector<std::pair<esglobal, esglobal> > >(mesh.neighbours().size()));
//	std::vector<std::vector<std::pair<esglobal, esglobal> > > rBuffer(mesh.neighbours().size());
//	std::vector<std::vector<std::pair<esglobal, esglobal> > > nepairs(threads);
//
//	#pragma omp parallel for
//	for (size_t t = 0; t < threads; t++) {
//		const crsiterator<eslocal>& enodes = pstore->nodesIndices->iterator(t);
//		const crsiterator<eslocal> nranks = chstore->ranks->iterator();
//		while (enodes.next()) {
//			for (auto n = enodes.begin(); n != enodes.end(); ++n) {
//				nepairs[t].push_back(std::make_pair(mesh.coordinates().globalIndex(*n), enodes.index() + pstore->elementOffset));
//
//				nranks.seek(*n);
//				for (auto c = nranks.begin(); c != nranks.end(); ++c) {
//					if (*c != environment->MPIrank) {
//						sBuffer[t][n2i(*c)].push_back(std::make_pair(mesh.coordinates().globalIndex(*n), enodes.index() + pstore->elementOffset));
//					}
//				}
//			}
//		}
//		std::sort(nepairs[t].begin(), nepairs[t].end());
//	}
//
//	#pragma omp parallel for
//	for (size_t n = 0; n < sBuffer[0].size(); n++) {
//		for (size_t t = 1; t < threads; t++) {
//			sBuffer[0][n].insert(sBuffer[0][n].end(), sBuffer[t][n].begin(), sBuffer[t][n].end());
//		}
//		std::sort(sBuffer[0][n].begin(), sBuffer[0][n].end());
//	}
//
//	if (!Communication::exchangeUnknownSize(sBuffer[0], rBuffer, mesh.neighbours())) {
//		ESINFO(ERROR) << "ESPRESO internal error while exchanging parent elements.";
//	}
//
//	std::vector<std::vector<std::pair<esglobal, esglobal> > > nparents(threads);
//
//	#pragma omp parallel for
//	for (size_t t = 0; t < threads; t++) {
//		if (chstore->distribution[t] == chstore->distribution[t + 1]) {
//			continue;
//		}
//		std::pair<esglobal, esglobal> first(mesh.coordinates().globalIndex(chstore->distribution[t]), 0);
//		std::pair<esglobal, esglobal>  last(mesh.coordinates().globalIndex(chstore->distribution[t + 1] - 1), 0);
//
//		auto compare = [] (const std::pair<esglobal, esglobal> &p, const std::pair<esglobal, esglobal> &node) { return p.first < node.first; };
//
//		for (size_t n = 0; n < rBuffer.size(); n++) {
//			auto begin = std::lower_bound(rBuffer[n].begin(), rBuffer[n].end(), first, compare);
//			auto end   = std::upper_bound(rBuffer[n].begin(), rBuffer[n].end(), last , compare);
//			if (begin < end) {
//				nparents[t].insert(nparents[t].end(), begin, end);
//			}
//		}
//		for (size_t p = 0; p < threads; p++) {
//			auto begin = std::lower_bound(nepairs[p].begin(), nepairs[p].end(), first, compare);
//			auto end   = std::upper_bound(nepairs[p].begin(), nepairs[p].end(), last , compare);
//			if (begin < end) {
//				nparents[t].insert(nparents[t].end(), begin, end);
//			}
//		}
//		std::sort(nparents[t].begin(), nparents[t].end());
//	}
//
//	std::vector<size_t> odistribution = { 0 };
//	std::vector<size_t> edistribution = { 0 };
//	for (size_t t = 0; t < threads; t++) {
//		odistribution.push_back(odistribution.back() + chstore->distribution[t + 1] - chstore->distribution[t]);
//		edistribution.push_back(edistribution.back() + nparents[t].size());
//	}
//	odistribution.back() += 1;
//
//	parray<esglobal> offsets(odistribution.size(), odistribution.data());
//	parray<esglobal> nelements(edistribution.size(), odistribution.data());
//
//	#pragma omp parallel for
//	for (size_t t = 0; t < threads; t++) {
//		size_t index = chstore->distribution[t];
//		size_t offset = edistribution[t];
//
//		offsets[index] = offset;
//		for (size_t n = 0; n < nparents[t].size(); n++, offset++) {
//			nelements[offset] = nparents[t][n].second;
//			if (n && nparents[t][n - 1].first != nparents[t][n].first) {
//				offsets[++index] = offset;
//			}
//		}
//	}
//	offsets[offsets.size() - 1] = edistribution.back();
//
//	chstore->parentIndices = new pcrsdata<esglobal>(std::move(offsets), std::move(nelements));
}

void Transformation::createFullDualGraph(NewMesh &mesh)
{
//	if (mesh.__nodes->parentIndices == NULL) {
//		Transformation::fillParentElements(mesh, TFlags::ELEVEL::NODE, TFlags::ELEVEL::ELEMENT);
//	}
//	if (mesh.__elements->haloElements == NULL && environment->MPIsize) {
//		Transformation::fillHaloElement(mesh);
//	}
//
//	size_t threads = environment->OMP_NUM_THREADS;
//
//	std::vector<std::vector<esglobal> > ensize(threads);
//	std::vector<std::vector<esglobal> > eneighbors(threads);
//
//	const parray<Element*> &element = *mesh.__elements->elements;
//	size_t eoffset = mesh.__elements->elementOffset;
//
//	#pragma omp parallel for
//	for (size_t t = 0; t < threads; t++) {
//		const crsiterator<esglobal> &eindices = mesh.__elements->nodesIndices->iterator(t);
//		const crsiterator<esglobal> nelements = mesh.__nodes->parentIndices->iterator();
//		std::vector<esglobal> neighbors;
//		size_t myCommon, otherCommon, first, ncounter;
//		HaloElement h;
//
//		while (eindices.next()) {
//			neighbors.clear();
//			for (auto n = eindices.begin(); n != eindices.end(); ++n) {
//				nelements.seek(*n);
//				neighbors.insert(neighbors.end(), nelements.begin(), nelements.end());
//			}
//			std::sort(neighbors.begin(), neighbors.end());
//			myCommon = element[eindices.index()]->nCommon();
//			first = 0;
//			ensize[t].push_back(0);
//			ncounter = 1;
//			while (first < neighbors.size()) {
//				while (first + ncounter < neighbors.size() && neighbors[first] == neighbors[first + ncounter]) {
//					ncounter++;
//				}
//				if (neighbors[first] != eindices.index() + eoffset) {
//					if (neighbors[first] - eoffset < mesh.__elements->elements->size()) {
//						otherCommon = element[neighbors[first] - eoffset]->nCommon();
//					} else {
//						h.id = neighbors[first];
//						otherCommon = std::lower_bound(mesh.__elements->haloElements->begin(), mesh.__elements->haloElements->end(), h)->nCommon;
//					}
//					if (ncounter >= std::min(myCommon, otherCommon)) {
//						eneighbors[t].push_back(neighbors[first]);
//						++ensize[t].back();
//					}
//				}
//				first += ncounter;
//				ncounter = 1;
//			}
//		}
//	}
//
//	mesh.__elements->fullDual = new pcrsdata<esglobal>(ensize, eneighbors);
}

void Transformation::createRestrictedDualGraph(NewMesh &mesh, TFlags::SEPARATE separate)
{
//	if (mesh.__nodes->parentIndices == NULL) {
//		Transformation::fillParentElements(mesh, TFlags::ELEVEL::NODE, TFlags::ELEVEL::ELEMENT);
//	}
//
//	size_t threads = environment->OMP_NUM_THREADS;
//
//	std::vector<std::vector<esglobal> > ensize(threads);
//	std::vector<std::vector<esglobal> > eneighbors(threads);
//
//	const parray<Element*> &element = *mesh.__elements->elements;
//	size_t eoffset = mesh.__elements->elementOffset;
//
//	#pragma omp parallel for
//	for (size_t t = 0; t < threads; t++) {
//		const crsiterator<esglobal> &eindices = mesh.__elements->nodesIndices->iterator(t);
//		const crsiterator<esglobal> parents = mesh.__nodes->parentIndices->iterator();
//		std::vector<esglobal> neighbors;
//		int body1, body2, etype1, etype2, material1, material2, common1, common2, first, ncounter;
//
//		while (eindices.next()) {
//			neighbors.clear();
//			for (auto n = eindices.begin(); n != eindices.end(); ++n) {
//				parents.seek(*n);
//				neighbors.insert(neighbors.end(), parents.begin(), parents.end());
//			}
//			std::sort(neighbors.begin(), neighbors.end());
//
//			common1 = element[eindices.index()]->nCommon();
//			body1 = element[eindices.index()]->param(Element::Params::BODY);
//			etype1 = (int)element[eindices.index()]->type();
//			material1 = element[eindices.index()]->param(Element::Params::MATERIAL);
//
//			first = 0;
//			ensize[t].push_back(0);
//			ncounter = 1;
//			while (first < neighbors.size()) {
//				while (first + ncounter < neighbors.size() && neighbors[first] == neighbors[first + ncounter]) {
//					ncounter++;
//				}
//				if (neighbors[first] != eindices.index() + eoffset) {
//					if (neighbors[first] - eoffset < mesh.__elements->elements->size()) {
//						common2 = element[neighbors[first] - eoffset]->nCommon();
//						body2 = element[neighbors[first] - eoffset]->param(Element::Params::BODY);
//						etype2 = (int)element[neighbors[first] - eoffset]->type();
//						material2 = element[neighbors[first] - eoffset]->param(Element::Params::MATERIAL);
//
//						if (
//								(ncounter >= std::min(common1, common2)) &&
//								(body1 == body2) &&
//								(etype1 == etype2 || !(separate & TFlags::SEPARATE::DEGREES_OF_FREEDOM)) &&
//								(material1 == material2 || !(separate & TFlags::SEPARATE::MATERIALS))) {
//
//							eneighbors[t].push_back(neighbors[first] - eoffset);
//							++ensize[t].back();
//						}
//					}
//				}
//				first += ncounter;
//				ncounter = 1;
//			}
//		}
//	}
//
//	mesh.__elements->restrictedDual = new pcrsdata<esglobal>(ensize, eneighbors);
}



