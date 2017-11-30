
#ifndef SRC_MESH_ELEMENTS_PLANE_SQUARE4_H_
#define SRC_MESH_ELEMENTS_PLANE_SQUARE4_H_

#include "../element.h"

namespace espreso {

struct Square4 {

	static Element fill(Element e, Element* begin);

	static Element create(Element* begin)
	{
		return fill(Element(Element::TYPE::PLANE, Element::CODE::SQUARE4, 4, 2, 1), begin);
	}
};

}


#endif /* SRC_MESH_ELEMENTS_PLANE_SQUARE4_H_ */
