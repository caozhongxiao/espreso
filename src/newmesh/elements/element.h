
#ifndef SRC_NEWMESH_ELEMENTS_NEWELEMENT_H_
#define SRC_NEWMESH_ELEMENTS_NEWELEMENT_H_

#include <cstddef>
#include <vector>

namespace espreso {

template <typename TEBoundaries, typename TEData> class serializededata;
class ElementStore;

struct Element {

	enum class TYPE: int {
		POINT  = 0,
		LINE   = 1,
		PLANE  = 2,
		VOLUME = 3,
	};

	enum class CODE: int {
		POINT1,

		// without mid-points
		LINE2,

		TRIANGLE3,
		SQUARE4,

		TETRA4,
		PYRAMID5,
		PRISMA6,
		HEXA8,

		// with mid-points
		LINE3,

		TRIANGLE6,
		SQUARE8,

		TETRA10,
		PYRAMID13,
		PRISMA15,
		HEXA20,

		// number of element types
		SIZE
	};

	TYPE type;
	CODE code;
	int coarseNodes;
	int nCommonFace;
	int nCommonEdge;
	// add base functions

	serializededata<int, int> *faces;
	serializededata<int, int> *edges;

	serializededata<int, Element*> *facepointers;
	serializededata<int, Element*> *edgepointers;

	Element(TYPE type, CODE code, int coarseNodes, int nCommonFace, int nCommonEdge)
	: type(type), code(code), coarseNodes(coarseNodes), nCommonFace(nCommonFace), nCommonEdge(nCommonEdge),
	  faces(NULL), edges(NULL), facepointers(NULL), edgepointers(NULL) {}

	Element()
	: type(TYPE::POINT), code(CODE::POINT1), coarseNodes(1), nCommonFace(1), nCommonEdge(1),
	  faces(NULL), edges(NULL), facepointers(NULL), edgepointers(NULL) {}
};

}



#endif /* SRC_NEWMESH_ELEMENTS_NEWELEMENT_H_ */
