
#ifndef SRC_MESH_ELEMENTS_PLANE_SQUARE8_H_
#define SRC_MESH_ELEMENTS_PLANE_SQUARE8_H_

#include "mesh/elements/element.h"

namespace espreso {

struct Square8 {

	static Element fill(Element e, Element* begin);

	static Element create(Element* begin)
	{
		return fill(Element(Element::TYPE::PLANE, Element::CODE::SQUARE8, 8, 4, 3, 2), begin);
	}
};

}


#endif /* SRC_MESH_ELEMENTS_PLANE_SQUARE8_H_ */
