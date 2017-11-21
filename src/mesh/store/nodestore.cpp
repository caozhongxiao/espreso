
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
  ranks(NULL),

  idomains(NULL)
{

}

NodeStore::~NodeStore()
{
	if (IDs == NULL) { delete IDs; }
	if (elements == NULL) { delete elements; }

	if (coordinates == NULL) { delete coordinates; }
	if (ranks == NULL) { delete ranks; }

	if (idomains == NULL) { delete idomains; }
}

void NodeStore::store(const std::string &file)
{
	std::ofstream os(file + std::to_string(environment->MPIrank) + ".txt");

	Store::storedata(os, "IDs", IDs);
	Store::storedata(os, "elements", elements);

	Store::storedata(os, "coordinates", coordinates);
	Store::storedata(os, "ranks", ranks);

	Store::storedata(os, "ineighbors", idomains);
}

void NodeStore::permute(const std::vector<eslocal> &permutation, const std::vector<size_t> &distribution)
{
	this->distribution = distribution;

	if (IDs != NULL) { IDs->permute(permutation, distribution); }
	if (elements != NULL) { elements->permute(permutation, distribution); }

	if (coordinates != NULL) { coordinates->permute(permutation, distribution); }
	if (ranks != NULL) { ranks->permute(permutation, distribution); }
}

std::vector<eslocal> NodeStore::gatherNodeDistribution()
{
	return Store::gatherDistribution(size);
}



