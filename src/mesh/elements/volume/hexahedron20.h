
#ifndef SRC_MESH_ELEMENTS_VOLUME_HEXAHEDRON20_H_
#define SRC_MESH_ELEMENTS_VOLUME_HEXAHEDRON20_H_

#include "mesh/elements/element.h"

namespace espreso {

struct Hexahedron20 {

	static Element fill(Element e, Element* begin);

	static Element create(Element* begin)
	{
		return fill(Element(Element::TYPE::VOLUME, Element::CODE::HEXA20, 20, 8, 4, 3), begin);
	}
};

}


#endif /* SRC_MESH_ELEMENTS_VOLUME_HEXAHEDRON20_H_ */
