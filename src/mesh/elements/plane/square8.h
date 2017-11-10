
#ifndef SRC_MESH_ELEMENTS_PLANE_SQUARE8_H_
#define SRC_MESH_ELEMENTS_PLANE_SQUARE8_H_

#include "../element.h"

namespace espreso {

struct Square8 {

	static Element create()
	{
		return Element(Element::TYPE::PLANE, Element::CODE::SQUARE8, 8, 3, 2);
	}
};

}


#endif /* SRC_MESH_ELEMENTS_PLANE_SQUARE8_H_ */
