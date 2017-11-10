
#ifndef SRC_NEWMESH_ELEMENTS_PLANE_SQUARE4_H_
#define SRC_NEWMESH_ELEMENTS_PLANE_SQUARE4_H_

#include "../element.h"

namespace espreso {

struct Square4 {

	static Element create()
	{
		return Element(Element::TYPE::PLANE, Element::CODE::SQUARE4, 4, 2, 1);
	}
};

}


#endif /* SRC_NEWMESH_ELEMENTS_PLANE_SQUARE4_H_ */
