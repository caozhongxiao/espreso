
#ifndef SRC_MESH_ELEMENTS_VOLUME_TETRAHEDRON4_H_
#define SRC_MESH_ELEMENTS_VOLUME_TETRAHEDRON4_H_

#include "../element.h"

namespace espreso {

struct Tetrahedron4 {

	static Element fill(Element e, size_t thread, Element* begin);

	static Element create(size_t thread, Element* begin)
	{
		return fill(Element(Element::TYPE::VOLUME, Element::CODE::TETRA4, 4, 3, 2), thread, begin);
	}
};

}


#endif /* SRC_MESH_ELEMENTS_VOLUME_TETRAHEDRON4_H_ */
