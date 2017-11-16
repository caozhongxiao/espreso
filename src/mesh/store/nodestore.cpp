
#include "store.h"
#include "nodestore.h"

#include "../../basis/containers/point.h"
#include "../../basis/containers/serializededata.h"

using namespace espreso;

NodeStore::NodeStore()
: size(0),
  distribution({0, 0}),

  IDs(NULL),
  elements(NULL),

  coordinates(NULL),
  domains(NULL),
  ranks(NULL)
{

}

NodeStore::~NodeStore()
{
	if (IDs == NULL) { delete IDs; }
	if (elements == NULL) { delete elements; }

	if (coordinates == NULL) { delete coordinates; }
	if (domains == NULL) { delete domains; }
	if (ranks == NULL) { delete ranks; }
}

void NodeStore::store(const std::string &file)
{
	std::ofstream os(file + std::to_string(environment->MPIrank) + ".txt");

	Store::storedata(os, "IDs", IDs);
	Store::storedata(os, "elements", elements);

	Store::storedata(os, "coordinates", coordinates);
	Store::storedata(os, "domains", domains);
	Store::storedata(os, "ranks", ranks);
}

void NodeStore::permute(const std::vector<eslocal> &permutation, const std::vector<size_t> &distribution)
{
	this->distribution = distribution;

	if (IDs != NULL) { IDs->permute(permutation, distribution); }
	if (elements != NULL) { elements->permute(permutation, distribution); }

	if (coordinates != NULL) { coordinates->permute(permutation, distribution); }
	if (domains != NULL) { domains->permute(permutation, distribution); }
	if (ranks != NULL) { ranks->permute(permutation, distribution); }
}

std::vector<eslocal> NodeStore::gatherNodeDistribution()
{
	return Store::gatherDistribution(size);
}



