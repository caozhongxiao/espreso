
#ifndef SRC_MESH_ELEMENTS_PLANE_TRIANGLE3_H_
#define SRC_MESH_ELEMENTS_PLANE_TRIANGLE3_H_

#include "../element.h"

namespace espreso {

struct Triangle3 {

	static Element fill(Element e, Element* begin);

	static Element create(Element* begin)
	{
		return fill(Element(Element::TYPE::PLANE, Element::CODE::TRIANGLE3, 3, 2, 1), begin);
	}
};

}


#endif /* SRC_MESH_ELEMENTS_PLANE_TRIANGLE3_H_ */
