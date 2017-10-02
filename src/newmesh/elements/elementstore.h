
#ifndef SRC_NEWMESH_ELEMENTS_ELEMENTSTORE_H_
#define SRC_NEWMESH_ELEMENTS_ELEMENTSTORE_H_

namespace espreso {

template <typename TEBoundaries, typename TEData> class serializededata;
struct Point;

struct ElementStore {

	serializededata<eslocal, esglobal>* IDs;

	serializededata<eslocal, eslocal>* elems;
	serializededata<eslocal, eslocal>* faces;
	serializededata<eslocal, eslocal>* edges;
	serializededata<eslocal, eslocal>* nodes;

	serializededata<eslocal, Point>* coordinates;

	serializededata<eslocal, esglobal>* globalDual;
	serializededata<eslocal, eslocal>* localDual;

	ElementStore();
	~ElementStore();
};

}



#endif /* SRC_NEWMESH_ELEMENTS_ELEMENTSTORE_H_ */
