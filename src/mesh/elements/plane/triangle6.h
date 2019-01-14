
#ifndef SRC_MESH_ELEMENTS_PLANE_TRIANGLE6_H_
#define SRC_MESH_ELEMENTS_PLANE_TRIANGLE6_H_

#include "mesh/elements/element.h"

namespace espreso {

struct Triangle6 {

	static Element fill(Element e, Element* begin);

	static Element create(Element* begin)
	{
		return fill(Element(Element::TYPE::PLANE, Element::CODE::TRIANGLE6, 6, 3, 3, 2), begin);
	}
};

}


#endif /* SRC_MESH_ELEMENTS_PLANE_TRIANGLE6_H_ */
