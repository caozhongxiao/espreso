
#ifndef SRC_MESH_ELEMENTS_VOLUME_TETRAHEDRON4_H_
#define SRC_MESH_ELEMENTS_VOLUME_TETRAHEDRON4_H_

#include "../element.h"

namespace espreso {

struct Tetrahedron4 {

	static Element fill(Element e, Element* begin);

	static Element create(Element* begin)
	{
		return fill(Element(Element::TYPE::VOLUME, Element::CODE::TETRA4, 4, 4, 3, 2), begin);
	}
};

}


#endif /* SRC_MESH_ELEMENTS_VOLUME_TETRAHEDRON4_H_ */
