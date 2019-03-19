
#include "surfacestore.h"

#include "basis/containers/point.h"
#include "basis/containers/serializededata.h"
#include "basis/utilities/packing.h"
#include "esinfo/envinfo.h"

#include "mesh/mesh.h"
#include "mesh/elements/element.h"

using namespace espreso;

SurfaceStore::SurfaceStore()
: triangles(NULL),
  elements(NULL),
  nelements(NULL),
  IDs(NULL),
  neighbors(NULL),
  nodes(NULL),
  coordinates(NULL),
  eoffset(0),
  epointers(NULL)
{

}

SurfaceStore::~SurfaceStore()
{
	if (triangles != NULL && triangles != elements) { delete triangles; }
	if (elements != NULL) { delete elements; }
	if (nelements != NULL) { delete nelements; }
	if (IDs != NULL) { delete IDs; }
	if (neighbors != NULL) { delete neighbors; }
	if (nodes != NULL) { delete nodes; }
	if (coordinates != NULL) { delete coordinates; }
	if (epointers != NULL) { delete epointers; }
}

size_t SurfaceStore::packedFullSize() const
{
	size_t packedSize = 0;

	packedSize += utils::packedSize(triangles);
	packedSize += utils::packedSize(elements);
	packedSize += utils::packedSize(nelements);
	packedSize += utils::packedSize(IDs);
	packedSize += utils::packedSize(neighbors);
	packedSize += utils::packedSize(nodes);
	packedSize += utils::packedSize(coordinates);

	packedSize += utils::packedSize(tdistribution);
	packedSize += utils::packedSize(edistribution);
	packedSize += utils::packedSize(cdistribution);

	packedSize += utils::packedSize(eoffset);
	packedSize += 1;
	if (epointers != NULL) {
		packedSize += sizeof(size_t) + epointers->datatarray().size() * sizeof(int);
	}
	packedSize += utils::packedSize(ecounters);

	return packedSize;
}

void SurfaceStore::packFull(char* &p) const
{
	utils::pack(triangles, p);
	utils::pack(elements, p);
	utils::pack(nelements, p);
	utils::pack(IDs, p);
	utils::pack(neighbors, p);
	utils::pack(nodes, p);
	utils::pack(coordinates, p);

	utils::pack(tdistribution, p);
	utils::pack(edistribution, p);
	utils::pack(cdistribution, p);

	utils::pack(eoffset, p);
	utils::pack(epointers != NULL, p);
	if (epointers != NULL) {
		std::vector<int> eindices;
		eindices.reserve(epointers->datatarray().size());
		for (size_t i = 0; i < epointers->datatarray().size(); ++i) {
			eindices.push_back(epointers->datatarray()[i] - Mesh::edata);
		}
		utils::pack(eindices, p);
	}
	utils::pack(ecounters, p);
}

void SurfaceStore::unpackFull(const char* &p)
{
	utils::unpack(triangles, p);
	utils::unpack(elements, p);
	utils::unpack(nelements, p);
	utils::unpack(IDs, p);
	utils::unpack(neighbors, p);
	utils::unpack(nodes, p);
	utils::unpack(coordinates, p);

	utils::unpack(tdistribution, p);
	utils::unpack(edistribution, p);
	utils::unpack(cdistribution, p);

	utils::unpack(eoffset, p);
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
	utils::unpack(ecounters, p);
}

