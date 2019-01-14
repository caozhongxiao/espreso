
#ifndef SRC_MESH_ELEMENTS_VOLUME_TETRAHEDRON10_H_
#define SRC_MESH_ELEMENTS_VOLUME_TETRAHEDRON10_H_

#include "mesh/elements/element.h"

namespace espreso {

struct Tetrahedron10 {

	static Element fill(Element e, Element* begin);

	static Element create(Element* begin)
	{
		return fill(Element(Element::TYPE::VOLUME, Element::CODE::TETRA10, 10, 4, 4, 3), begin);
	}
};

}

#endif /* SRC_MESH_ELEMENTS_VOLUME_TETRAHEDRON10_H_ */
