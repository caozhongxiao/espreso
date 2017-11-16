
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

  _eclasses(eclasses)
{

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

std::vector<eslocal> ElementStore::gatherElementProcDistribution()
{
	return Store::gatherDistribution(size);
}
