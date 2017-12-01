
#include "boundaryregionstore.h"

#include "../../basis/containers/serializededata.h"
#include "../../basis/utilities/utils.h"

#include "../../config/ecf/environment.h"

#include "../../mesh/elements/element.h"

using namespace espreso;

BoundaryRegionStore::BoundaryRegionStore(const std::string &name, std::vector<Element*> &eclasses)
: name(name),

  uniqueOffset(0),
  uniqueSize(0),
  uniqueTotalSize(0),

  faces(NULL),
  edges(NULL),
  nodes(NULL),

  facepointers(NULL),
  edgepointers(NULL),

  _eclasses(eclasses)
{

}

BoundaryRegionStore::~BoundaryRegionStore()
{
	if (faces == NULL) { delete faces; }
	if (edges == NULL) { delete edges; }
	if (nodes == NULL) { delete nodes; }

	if (facepointers == NULL) { delete facepointers; }
	if (edgepointers == NULL) { delete edgepointers; }
}

size_t BoundaryRegionStore::packedSize() const
{
	size_t size = Esutils::packedSize(name) + Esutils::packedSize(uniqueOffset) + Esutils::packedSize(uniqueSize) + Esutils::packedSize(uniqueTotalSize);
	size += sizeof(int); // bit-mask of what we need to pack
	if (faces) {
		size += faces->packedSize();
		size += sizeof(size_t) + facepointers->datatarray().size() * sizeof(int);
		size += Esutils::packedSize(facesIntervals);
	}
	if (edges) {
		size += edges->packedSize();
		size += sizeof(size_t) + edgepointers->datatarray().size() * sizeof(int);
		size += Esutils::packedSize(edgesIntervals);
	}
	if (nodes) {
		size += nodes->packedSize();
		size += Esutils::packedSize(nodesIntervals);
	}
	return size;
}

void BoundaryRegionStore::pack(char* &p) const
{
	Esutils::pack(name, p);
	Esutils::pack(uniqueOffset, p);
	Esutils::pack(uniqueSize, p);
	Esutils::pack(uniqueTotalSize, p);
	int bitmask = 0;
	if (faces) { bitmask += 4; }
	if (edges) { bitmask += 2; }
	if (nodes) { bitmask += 1; }
	Esutils::pack(bitmask, p);
	if (faces) {
		faces->pack(p);
		std::vector<int> eindices;
		eindices.reserve(facepointers->datatarray().size());

		size_t threads = environment->OMP_NUM_THREADS;
		for (size_t t = 0; t < threads; t++) {
			for (size_t i = facepointers->datatarray().distribution()[t]; i < facepointers->datatarray().distribution()[t + 1]; ++i) {
				eindices.push_back(facepointers->datatarray()[i] - _eclasses[t]);
			}
		}
		Esutils::pack(eindices, p);
		Esutils::pack(facesIntervals, p);
	}
	if (edges) {
		edges->pack(p);
		std::vector<int> eindices;
		eindices.reserve(edgepointers->datatarray().size());

		size_t threads = environment->OMP_NUM_THREADS;
		for (size_t t = 0; t < threads; t++) {
			for (size_t i = edgepointers->datatarray().distribution()[t]; i < edgepointers->datatarray().distribution()[t + 1]; ++i) {
				eindices.push_back(edgepointers->datatarray()[i] - _eclasses[t]);
			}
		}
		Esutils::pack(eindices, p);
		Esutils::pack(edgesIntervals, p);
	}
	if (nodes) {
		nodes->pack(p);
		Esutils::pack(nodesIntervals, p);
	}
}

void BoundaryRegionStore::unpack(const char* &p)
{
	Esutils::unpack(name, p);
	Esutils::unpack(uniqueOffset, p);
	Esutils::unpack(uniqueSize, p);
	Esutils::unpack(uniqueTotalSize, p);
	int bitmask;
	Esutils::unpack(bitmask, p);

	if (bitmask & 4) {
		if (faces == NULL) {
			faces = new serializededata<eslocal, eslocal>(tarray<eslocal>(1, 0), tarray<eslocal>(1, 0));
			facepointers = new serializededata<eslocal, Element*>(1, tarray<Element*>(1, 0));
		}
		faces->unpack(p);
		std::vector<int> eindices;
		Esutils::unpack(eindices, p);
		if (facepointers != NULL) {
			delete facepointers;
		}
		facepointers = new serializededata<eslocal, Element*>(1, tarray<Element*>(1, eindices.size()));
		for (size_t i = 0; i < eindices.size(); ++i) {
			facepointers->datatarray()[i] = &_eclasses[0][eindices[i]];
		}
		Esutils::unpack(facesIntervals, p);
	}

	if (bitmask & 2) {
		if (edges == NULL) {
			edges = new serializededata<eslocal, eslocal>(tarray<eslocal>(1, 0), tarray<eslocal>(1, 0));
			edgepointers = new serializededata<eslocal, Element*>(1, tarray<Element*>(1, 0));
		}
		edges->unpack(p);
		std::vector<int> eindices;
		Esutils::unpack(eindices, p);
		if (edgepointers != NULL) {
			delete edgepointers;
		}
		edgepointers = new serializededata<eslocal, Element*>(1, tarray<Element*>(1, eindices.size()));
		for (size_t i = 0; i < eindices.size(); ++i) {
			edgepointers->datatarray()[i] = &_eclasses[0][eindices[i]];
		}
		Esutils::unpack(edgesIntervals, p);
	}

	if (bitmask & 1) {
		if (nodes == NULL) {
			nodes = new serializededata<eslocal, eslocal>(tarray<eslocal>(1, 0), tarray<eslocal>(1, 0));
		}
		nodes->unpack(p);
		Esutils::unpack(nodesIntervals, p);
	}
}
