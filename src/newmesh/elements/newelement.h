
#ifndef SRC_NEWMESH_ELEMENTS_NEWELEMENT_H_
#define SRC_NEWMESH_ELEMENTS_NEWELEMENT_H_

#include <cstddef>

namespace espreso {

template <typename TEData> class edata;
class ElementStore;

struct NewElement {

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
		HEXA20
	};

	TYPE type;
	CODE code;
	int coarseNodes;
	int nCommonFace;
	int nCommonEdge;
	// add base functions


	NewElement(TYPE type, CODE code, int coarseNodes, int nCommonFace, int nCommonEdge)
	: type(type), code(code), coarseNodes(coarseNodes), nCommonFace(nCommonFace), nCommonEdge(nCommonEdge) {}

	NewElement()
	: type(TYPE::POINT), code(CODE::POINT1), coarseNodes(1), nCommonFace(1), nCommonEdge(1) {}
};

}



#endif /* SRC_NEWMESH_ELEMENTS_NEWELEMENT_H_ */
