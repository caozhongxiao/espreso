
#ifndef SRC_NEWMESH_ELEMENTS_PLANE_TRIANGLE6_H_
#define SRC_NEWMESH_ELEMENTS_PLANE_TRIANGLE6_H_

#include "../element.h"

namespace espreso {

struct Triangle6 {

	static Element create()
	{
		return Element(Element::TYPE::PLANE, Element::CODE::TRIANGLE6, 6, 3, 2);
	}
};

}


#endif /* SRC_NEWMESH_ELEMENTS_PLANE_TRIANGLE6_H_ */
