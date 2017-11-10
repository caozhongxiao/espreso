
#ifndef SRC_NEWMESH_ELEMENTS_PLANE_TRIANGLE3_H_
#define SRC_NEWMESH_ELEMENTS_PLANE_TRIANGLE3_H_

#include "../element.h"

namespace espreso {

struct Triangle3 {

	static Element create()
	{
		return Element(Element::TYPE::PLANE, Element::CODE::TRIANGLE3, 3, 2, 1);
	}
};

}


#endif /* SRC_NEWMESH_ELEMENTS_PLANE_TRIANGLE3_H_ */
