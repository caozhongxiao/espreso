
#ifndef SRC_MESH_ELEMENTS_ELEMENT_H_
#define SRC_MESH_ELEMENTS_ELEMENT_H_

#include <vector>

namespace espreso {

class DenseMatrix;
template <typename TEBoundaries, typename TEData> class serializededata;

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

	std::vector<DenseMatrix> *N;
	std::vector<DenseMatrix> *dN;
	std::vector<double> *weighFactor;

	serializededata<int, int> *faces;
	serializededata<int, int> *edges;

	serializededata<int, Element*> *facepointers;
	serializededata<int, Element*> *edgepointers;

	Element(TYPE type, CODE code, int coarseNodes, int nCommonFace, int nCommonEdge)
	: type(type), code(code), coarseNodes(coarseNodes), nCommonFace(nCommonFace), nCommonEdge(nCommonEdge),
	  N(NULL), dN(NULL), weighFactor(NULL),
	  faces(NULL), edges(NULL), facepointers(NULL), edgepointers(NULL) {}

	Element()
	: type(TYPE::POINT), code(CODE::POINT1), coarseNodes(1), nCommonFace(1), nCommonEdge(1),
	  N(NULL), dN(NULL), weighFactor(NULL),
	  faces(NULL), edges(NULL), facepointers(NULL), edgepointers(NULL) {}
};

}



#endif /* SRC_MESH_ELEMENTS_ELEMENT_H_ */
