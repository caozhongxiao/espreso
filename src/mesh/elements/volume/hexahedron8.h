
#ifndef SRC_MESH_ELEMENTS_VOLUME_HEXAHEDRON8_H_
#define SRC_MESH_ELEMENTS_VOLUME_HEXAHEDRON8_H_

#include "../element.h"

namespace espreso {

struct Hexahedron8 {

	static Element fill(Element e, Element* begin);

	static Element create(Element* begin)
	{
		return fill(Element(Element::TYPE::VOLUME, Element::CODE::HEXA8, 8, 8, 3, 2), begin);
	}

};

}


#endif /* SRC_MESH_ELEMENTS_VOLUME_HEXAHEDRON8_H_ */
