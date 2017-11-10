
#ifndef SRC_NEWMESH_ELEMENTS_VOLUME_HEXAHEDRON8_H_
#define SRC_NEWMESH_ELEMENTS_VOLUME_HEXAHEDRON8_H_

#include "../element.h"

namespace espreso {

struct Hexahedron8 {

	static Element fill(Element e, size_t thread, Element* begin);

	static Element create(size_t thread, Element* begin)
	{
		return fill(Element(Element::TYPE::VOLUME, Element::CODE::HEXA8, 8, 3, 2), thread, begin);
	}

};

}


#endif /* SRC_NEWMESH_ELEMENTS_VOLUME_HEXAHEDRON8_H_ */
