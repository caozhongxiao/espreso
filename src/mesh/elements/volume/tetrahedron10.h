
#ifndef SRC_MESH_ELEMENTS_VOLUME_TETRAHEDRON10_H_
#define SRC_MESH_ELEMENTS_VOLUME_TETRAHEDRON10_H_

#include "../element.h"

namespace espreso {

struct Tetrahedron10 {

	static Element create()
	{
		return Element(Element::TYPE::VOLUME, Element::CODE::TETRA10, 10, 4, 3);
	}
};

}

#endif /* SRC_MESH_ELEMENTS_VOLUME_TETRAHEDRON10_H_ */
