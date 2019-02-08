
#include "boundaryregionstore.h"

#include "esinfo/envinfo.h"
#include "basis/containers/serializededata.h"
#include "basis/utilities/utils.h"

#include "mesh/mesh.h"
#include "mesh/elements/element.h"

using namespace espreso;

BoundaryRegionStore::BoundaryRegionStore(const std::string &name)
: name(name),
  dimension(0),
  area(0),

  uniqueOffset(0),
  uniqueSize(0),
  uniqueTotalSize(0),

  procNodes(NULL),
  triangles(NULL),
  nodes(NULL),

  epointers(NULL),
  ecounters(static_cast<int>(Element::CODE::SIZE))
{

}

BoundaryRegionStore::~BoundaryRegionStore()
{
	if (procNodes == NULL) { delete procNodes; }
	if (triangles != NULL && triangles != procNodes) { delete triangles; }
	if (nodes == NULL) { delete nodes; }

	if (epointers == NULL) { delete epointers; }
}

void BoundaryRegionStore::permute(const std::vector<esint> &permutation, const std::vector<size_t> &distribution)
{
	this->distribution = distribution;

	if (procNodes != NULL) {
		procNodes->permute(permutation, distribution);
	}

	if (epointers != NULL) {
		epointers->permute(permutation, distribution);
	}

	eintervals.clear();
}

size_t BoundaryRegionStore::packedSize() const
{
	return
			Esutils::packedSize(name) +
			Esutils::packedSize(dimension) +
			Esutils::packedSize(uniqueOffset) +
			Esutils::packedSize(uniqueSize) +
			Esutils::packedSize(uniqueTotalSize) +
			procNodes->packedSize() +
			sizeof(size_t) + epointers->datatarray().size() * sizeof(int) +
			Esutils::packedSize(ecounters) +
			Esutils::packedSize(eintervals) +
			nodes->packedSize() +
			Esutils::packedSize(nintervals);
}

void BoundaryRegionStore::pack(char* &p) const
{
	Esutils::pack(name, p);
	Esutils::pack(dimension, p);
	Esutils::pack(uniqueOffset, p);
	Esutils::pack(uniqueSize, p);
	Esutils::pack(uniqueTotalSize, p);
	procNodes->pack(p);
	std::vector<int> eindices;
	eindices.reserve(epointers->datatarray().size());

	size_t threads = info::env::OMP_NUM_THREADS;
	for (size_t t = 0; t < threads; t++) {
		for (size_t i = epointers->datatarray().distribution()[t]; i < epointers->datatarray().distribution()[t + 1]; ++i) {
			eindices.push_back(epointers->datatarray()[i] - Mesh::edata);
		}
	}
	Esutils::pack(eindices, p);
	Esutils::pack(ecounters, p);
	Esutils::pack(eintervals, p);
	nodes->pack(p);
	Esutils::pack(nintervals, p);
}

void BoundaryRegionStore::unpack(const char* &p)
{
	Esutils::unpack(name, p);
	Esutils::unpack(dimension, p);
	Esutils::unpack(uniqueOffset, p);
	Esutils::unpack(uniqueSize, p);
	Esutils::unpack(uniqueTotalSize, p);
	if (procNodes == NULL) {
		procNodes = new serializededata<esint, esint>(tarray<esint>(1, 0), tarray<esint>(1, 0));
		epointers = new serializededata<esint, Element*>(1, tarray<Element*>(1, 0));
	}
	procNodes->unpack(p);
	std::vector<int> eindices;
	Esutils::unpack(eindices, p);
	if (epointers != NULL) {
		delete epointers;
	}
	epointers = new serializededata<esint, Element*>(1, tarray<Element*>(1, eindices.size()));
	for (size_t i = 0; i < eindices.size(); ++i) {
		epointers->datatarray()[i] = Mesh::edata + eindices[i];
	}
	Esutils::unpack(ecounters, p);
	Esutils::unpack(eintervals, p);
	if (nodes == NULL) {
		nodes = new serializededata<esint, esint>(tarray<esint>(1, 0), tarray<esint>(1, 0));
	}
	nodes->unpack(p);
	Esutils::unpack(nintervals, p);
}
