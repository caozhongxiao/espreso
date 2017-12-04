
#include "store.h"
#include "elementstore.h"

#include "../elements/element.h"

#include "../../basis/containers/serializededata.h"
#include "../../config/ecf/environment.h"

using namespace espreso;

ElementStore::ElementStore(std::vector<Element*> &eclasses)
: size(0),
  distribution({0, 0}),

  IDs(NULL),
  nodes(NULL),

  body(NULL),
  material(NULL),
  epointers(NULL),

  dual(NULL),
  decomposedDual(NULL),

  firstDomain(0),
  ndomains(1),
  nclusters(1),

  ecounters(static_cast<int>(Element::CODE::SIZE)),

  _eclasses(eclasses)
{

}

size_t ElementStore::packedSize() const
{
	if (nodes == NULL || epointers == NULL) {
		ESINFO(ERROR) << "ESPRESO internal error: invalid request for packedSize.";
	}
	return
			Esutils::packedSize(size) +
			nodes->packedSize() +
			sizeof(size_t) + epointers->datatarray().size() * sizeof(int) +
			Esutils::packedSize(firstDomain) +
			Esutils::packedSize(ndomains) +
			Esutils::packedSize(nclusters) +
			Esutils::packedSize(clusters) +
			Esutils::packedSize(elementsDistribution) +
			Esutils::packedSize(ecounters) +
			Esutils::packedSize(eintervals);
}

void ElementStore::pack(char* &p) const
{
	Esutils::pack(size, p);
	nodes->pack(p);
	if (epointers != NULL) {
		std::vector<int> eindices;
		eindices.reserve(epointers->datatarray().size());

		size_t threads = environment->OMP_NUM_THREADS;
		for (size_t t = 0; t < threads; t++) {
			for (size_t i = this->distribution[t]; i < this->distribution[t + 1]; ++i) {
				eindices.push_back(epointers->datatarray()[i] - _eclasses[t]);
			}
		}
		Esutils::pack(eindices, p);
	}
	Esutils::pack(firstDomain, p);
	Esutils::pack(ndomains, p);
	Esutils::pack(nclusters, p);
	Esutils::pack(clusters, p);
	Esutils::pack(elementsDistribution, p);
	Esutils::pack(ecounters, p);
	Esutils::pack(eintervals, p);
}

void ElementStore::unpack(const char* &p)
{
	if (nodes == NULL) {
		nodes = new serializededata<eslocal, eslocal>(tarray<eslocal>(1, 0), tarray<eslocal>(1, 0));
	}
	if (epointers == NULL) {
		epointers = new serializededata<eslocal, Element*>(1, tarray<Element*>(1, 0));
	}

	Esutils::unpack(size, p);
	nodes->unpack(p);
	if (epointers != NULL) {
		std::vector<int> eindices;
		Esutils::unpack(eindices, p);
		if (epointers != NULL) {
			delete epointers;
		}
		epointers = new serializededata<eslocal, Element*>(1, tarray<Element*>(1, size));
		for (size_t i = 0; i < size; ++i) {
			epointers->datatarray()[i] = &_eclasses[0][eindices[i]];
		}
	}
	Esutils::unpack(firstDomain, p);
	Esutils::unpack(ndomains, p);
	Esutils::unpack(nclusters, p);
	Esutils::unpack(clusters, p);
	Esutils::unpack(elementsDistribution, p);
	Esutils::unpack(ecounters, p);
	Esutils::unpack(eintervals, p);
}

ElementStore::~ElementStore()
{
	if (IDs == NULL) { delete IDs; }
	if (nodes == NULL) { delete nodes; }

	if (body == NULL) { delete body; }
	if (material == NULL) { delete material; }
	if (epointers == NULL) { delete epointers; }

	if (dual == NULL) { delete dual; }
	if (decomposedDual == NULL) { delete decomposedDual; }
}

void ElementStore::store(const std::string &file)
{
	std::ofstream os(file + std::to_string(environment->MPIrank) + ".txt");

	Store::storedata(os, "IDs", IDs);

	Store::storedata(os, "nodes", nodes);

	Store::storedata(os, "body", body);
	Store::storedata(os, "material", material);
	Store::storedata(os, "epointers", epointers);

	Store::storedata(os, "dual", dual);
	Store::storedata(os, "decomposedDual", decomposedDual);
}

void ElementStore::permute(const std::vector<eslocal> &permutation, const std::vector<size_t> &distribution)
{
	this->distribution = distribution;

	if (IDs != NULL) { IDs->permute(permutation, distribution); }

	if (nodes != NULL) { nodes->permute(permutation, distribution); }

	if (body != NULL) { body->permute(permutation, distribution); }
	if (material != NULL) { material->permute(permutation, distribution); }

	if (epointers != NULL) {
		size_t threads = environment->OMP_NUM_THREADS;
		if (threads > 1) {
			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				for (size_t i = epointers->datatarray().distribution()[t]; i < epointers->datatarray().distribution()[t + 1]; ++i) {
					epointers->datatarray()[i] = _eclasses[0] + (epointers->datatarray()[i] - _eclasses[t]);
				}
			}

			epointers->permute(permutation, distribution);

			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				for (size_t i = this->distribution[t]; i < this->distribution[t + 1]; ++i) {
					epointers->datatarray()[i] = _eclasses[t] + (epointers->datatarray()[i] - _eclasses[0]);
				}
			}
		} else {
			epointers->permute(permutation, distribution);
		}
	}

	if (dual != NULL) { dual->permute(permutation, distribution); }
	if (decomposedDual != NULL) { decomposedDual->permute(permutation, distribution); }
}

std::vector<eslocal> ElementStore::gatherElementsProcDistribution()
{
	return Store::gatherDistribution(size);
}

std::vector<eslocal> ElementStore::gatherDomainsProcDistribution()
{
	return Store::gatherDistribution(ndomains);
}

std::vector<eslocal> ElementStore::gatherDomainsDistribution()
{
	return Store::gatherDistribution(domainDistribution, firstDomain);
}

std::vector<eslocal> ElementStore::gatherElementsDistribution()
{
	return Store::gatherDistribution(elementsDistribution, IDs->datatarray().front());
}

std::vector<eslocal> ElementStore::gatherClustersDistribution()
{
	return Store::gatherDistribution(nclusters);
}
