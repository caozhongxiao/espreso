
#include "boundaryregionstore.h"

#include "esinfo/envinfo.h"
#include "basis/containers/serializededata.h"
#include "basis/utilities/packing.h"

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

BoundaryRegionStore::BoundaryRegionStore(const char* &packedData)
: BoundaryRegionStore("")
{
	unpackFull(packedData);
}

BoundaryRegionStore::~BoundaryRegionStore()
{
	if (procNodes != NULL) { delete procNodes; }
	if (triangles != NULL && triangles != procNodes) { delete triangles; }
	if (nodes != NULL) { delete nodes; }

	if (epointers != NULL) { delete epointers; }
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

size_t BoundaryRegionStore::packedFullSize() const
{
	size_t packedSize = 0;

	packedSize += utils::packedSize(name);
	packedSize += utils::packedSize(distribution);

	packedSize += utils::packedSize(dimension);
	packedSize += utils::packedSize(area);
	packedSize += utils::packedSize(uniqueOffset);
	packedSize += utils::packedSize(uniqueSize);
	packedSize += utils::packedSize(uniqueTotalSize);

	packedSize += utils::packedSize(procNodes);
	packedSize += utils::packedSize(triangles);
	packedSize += utils::packedSize(nodes);

	packedSize += 1;
	if (epointers != NULL) {
		packedSize += sizeof(size_t) + epointers->datatarray().size() * sizeof(int);
	}

	packedSize += utils::packedSize(eintervals);
	packedSize += utils::packedSize(eintervalsDistribution);
	packedSize += utils::packedSize(nintervals);

	packedSize += utils::packedSize(ecounters);

	return packedSize;
}

void BoundaryRegionStore::packFull(char* &p) const
{
	utils::pack(name, p);
	utils::pack(distribution, p);

	utils::pack(dimension, p);
	utils::pack(area, p);
	utils::pack(uniqueOffset, p);
	utils::pack(uniqueSize, p);
	utils::pack(uniqueTotalSize, p);

	utils::pack(procNodes, p);
	utils::pack(triangles, p);
	utils::pack(nodes, p);

	utils::pack(epointers != NULL, p);
	if (epointers != NULL) {
		std::vector<int> eindices;
		eindices.reserve(epointers->datatarray().size());
		for (size_t i = 0; i < epointers->datatarray().size(); ++i) {
			eindices.push_back(epointers->datatarray()[i] - Mesh::edata);
		}
		utils::pack(eindices, p);
	}

	utils::pack(eintervals, p);
	utils::pack(eintervalsDistribution, p);
	utils::pack(nintervals, p);

	utils::pack(ecounters, p);
}

void BoundaryRegionStore::unpackFull(const char* &p)
{
	utils::unpack(name, p);
	utils::unpack(distribution, p);

	utils::unpack(dimension, p);
	utils::unpack(area, p);
	utils::unpack(uniqueOffset, p);
	utils::unpack(uniqueSize, p);
	utils::unpack(uniqueTotalSize, p);

	utils::unpack(procNodes, p);
	utils::unpack(triangles, p);
	utils::unpack(nodes, p);

	bool notnull;
	utils::unpack(notnull, p);
	if (notnull) {
		std::vector<int> eindices;
		utils::unpack(eindices, p);
		if (epointers != NULL) {
			delete epointers;
		}
		epointers = new serializededata<esint, Element*>(1, tarray<Element*>(info::env::OMP_NUM_THREADS, eindices.size()));
		for (size_t i = 0; i < eindices.size(); ++i) {
			epointers->datatarray()[i] = Mesh::edata + eindices[i];
		}
	}

	utils::unpack(eintervals, p);
	utils::unpack(eintervalsDistribution, p);
	utils::unpack(nintervals, p);

	utils::unpack(ecounters, p);
}

size_t BoundaryRegionStore::packedSize() const
{
	return
			utils::packedSize(name) +
			utils::packedSize(dimension) +
			utils::packedSize(uniqueOffset) +
			utils::packedSize(uniqueSize) +
			utils::packedSize(uniqueTotalSize) +
			procNodes->packedSize() +
			sizeof(size_t) + epointers->datatarray().size() * sizeof(int) +
			utils::packedSize(ecounters) +
			utils::packedSize(eintervals) +
			nodes->packedSize() +
			utils::packedSize(nintervals);
}

void BoundaryRegionStore::pack(char* &p) const
{
	utils::pack(name, p);
	utils::pack(dimension, p);
	utils::pack(uniqueOffset, p);
	utils::pack(uniqueSize, p);
	utils::pack(uniqueTotalSize, p);
	procNodes->pack(p);
	std::vector<int> eindices;
	eindices.reserve(epointers->datatarray().size());

	size_t threads = info::env::OMP_NUM_THREADS;
	for (size_t t = 0; t < threads; t++) {
		for (size_t i = epointers->datatarray().distribution()[t]; i < epointers->datatarray().distribution()[t + 1]; ++i) {
			eindices.push_back(epointers->datatarray()[i] - Mesh::edata);
		}
	}
	utils::pack(eindices, p);
	utils::pack(ecounters, p);
	utils::pack(eintervals, p);
	nodes->pack(p);
	utils::pack(nintervals, p);
}

void BoundaryRegionStore::unpack(const char* &p)
{
	utils::unpack(name, p);
	utils::unpack(dimension, p);
	utils::unpack(uniqueOffset, p);
	utils::unpack(uniqueSize, p);
	utils::unpack(uniqueTotalSize, p);
	if (procNodes == NULL) {
		procNodes = new serializededata<esint, esint>(tarray<esint>(1, 0), tarray<esint>(1, 0));
		epointers = new serializededata<esint, Element*>(1, tarray<Element*>(1, 0));
	}
	procNodes->unpack(p);
	std::vector<int> eindices;
	utils::unpack(eindices, p);
	if (epointers != NULL) {
		delete epointers;
	}
	epointers = new serializededata<esint, Element*>(1, tarray<Element*>(1, eindices.size()));
	for (size_t i = 0; i < eindices.size(); ++i) {
		epointers->datatarray()[i] = Mesh::edata + eindices[i];
	}
	utils::unpack(ecounters, p);
	utils::unpack(eintervals, p);
	if (nodes == NULL) {
		nodes = new serializededata<esint, esint>(tarray<esint>(1, 0), tarray<esint>(1, 0));
	}
	nodes->unpack(p);
	utils::unpack(nintervals, p);
}
