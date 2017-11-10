
#include "elementstore.h"
#include "element.h"
#include "../../config/ecf/environment.h"
#include "../../basis/containers/serializededata.h"
#include "../../basis/point/point.h"
#include "../../basis/logging/logging.h"
#include "../../basis/utilities/communication.h"
#include "../../basis/utilities/utils.h"

#include <fstream>
#include <algorithm>
#include <numeric>

using namespace espreso;

ElementStore::ElementStore(std::vector<Element*> &eclasses)
: size(0),
  distribution({0, 0}),
  IDs(NULL),

  elems(NULL),
  faces(NULL),
  edges(NULL),
  nodes(NULL),

  coordinates(NULL),
  body(NULL),
  material(NULL),
  epointers(NULL),
  domains(NULL),
  ranks(NULL),

  dual(NULL),
  decomposedDual(NULL),

  _eclasses(eclasses)
{

}

ElementStore::~ElementStore()
{
	if (IDs == NULL) { delete IDs; }

	if (elems == NULL) { delete elems; }
	if (faces == NULL) { delete faces; }
	if (edges == NULL) { delete edges; }
	if (nodes == NULL) { delete nodes; }

	if (coordinates == NULL) { delete coordinates; }
	if (body == NULL) { delete body; }
	if (material == NULL) { delete material; }
	if (epointers == NULL) { delete epointers; }
	if (domains == NULL) { delete domains; }
	if (ranks == NULL) { delete ranks; }

	if (dual == NULL) { delete dual; }
	if (decomposedDual == NULL) { delete decomposedDual; }
}

template <typename TBoundaries, typename TData>
static void storedata(std::ofstream &os, const std::string &head, const serializededata<TBoundaries, TData> *data)
{
	if (data == NULL) {
		return;
	}

	os << head << "\n";
	for (auto elem = data->begin(); elem != data->end(); ++elem) {
		os << "[ ";
		for (auto i = elem->begin(); i != elem->end(); ++i) {
			os << *i << " ";
		}
		os << "] ";
	}
	os << "\n";
}

void ElementStore::store(const std::string &file)
{
	std::ofstream os(file + std::to_string(environment->MPIrank) + ".txt");

	storedata(os, "IDs", IDs);

	storedata(os, "elements", elems);
	storedata(os, "faces", faces);
	storedata(os, "edges", edges);
	storedata(os, "nodes", nodes);

	storedata(os, "coordinates", coordinates);
	storedata(os, "body", body);
	storedata(os, "material", material);
	storedata(os, "domains", domains);
	storedata(os, "ranks", ranks);

	storedata(os, "dual", dual);
	storedata(os, "decomposedDual", decomposedDual);

	if (epointers != NULL) {
		os << "epointers\n";
		for (auto elem = epointers->begin(); elem != epointers->end(); ++elem) {
			os << "[ " << static_cast<int>(elem->front()->code) << "] ";
		}
		os << "\n";
	}
}

void ElementStore::sort()
{
	if (IDs == NULL) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: sort element store that has no IDs.";
	}
	std::vector<eslocal> permutation(IDs->datatarray().size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) { return IDs->datatarray().data()[i] < IDs->datatarray().data()[j]; });

	permute(permutation);
}

void ElementStore::permute(const std::vector<eslocal> &permutation, const std::vector<size_t> *distribution)
{
	if (distribution != NULL) {
		this->distribution = *distribution;
	}

	if (IDs != NULL) { IDs->permute(permutation, distribution); }

	if (elems != NULL) { elems->permute(permutation, distribution); }
	if (faces != NULL) { faces->permute(permutation, distribution); }
	if (edges != NULL) { edges->permute(permutation, distribution); }
	if (nodes != NULL) { nodes->permute(permutation, distribution); }

	if (coordinates != NULL) { coordinates->permute(permutation, distribution); }
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

	if (domains != NULL) { domains->permute(permutation, distribution); }
	if (ranks != NULL) { ranks->permute(permutation, distribution); }

	if (dual != NULL) { dual->permute(permutation, distribution); }
	if (decomposedDual != NULL) { decomposedDual->permute(permutation, distribution); }
}

std::vector<esglobal> ElementStore::gatherElementDistrubution()
{
	std::vector<esglobal> result(environment->MPIsize + 1);
	esglobal esize = size;
	Communication::exscan(esize);

	MPI_Allgather(&esize, sizeof(esglobal), MPI_BYTE, result.data(), sizeof(esglobal), MPI_BYTE, environment->MPICommunicator);
	result.back() = esize + size;
	MPI_Bcast(&result.back(), sizeof(esglobal), MPI_BYTE, environment->MPIsize - 1, environment->MPICommunicator);

	return result;
}
