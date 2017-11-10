
#ifndef SRC_MESH_ELEMENTS_VOLUME_HEXAHEDRON20_H_
#define SRC_MESH_ELEMENTS_VOLUME_HEXAHEDRON20_H_

#include "../element.h"

namespace espreso {

struct Hexahedron20 {

	static Element create()
	{
		return Element(Element::TYPE::VOLUME, Element::CODE::HEXA20, 20, 4, 3);
	}
};

}


#endif /* SRC_MESH_ELEMENTS_VOLUME_HEXAHEDRON20_H_ */
