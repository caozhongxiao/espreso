
#ifndef SRC_MESH_ELEMENTS_VOLUME_PYRAMID13_H_
#define SRC_MESH_ELEMENTS_VOLUME_PYRAMID13_H_

#include "../element.h"

namespace espreso {

struct Pyramid13 {

	static Element fill(Element e, Element* begin);

	static Element create(Element* begin)
	{
		return fill(Element(Element::TYPE::VOLUME, Element::CODE::PYRAMID13, 13, 4, 3), begin);
	}
};

}


#endif /* SRC_MESH_ELEMENTS_VOLUME_PYRAMID13_H_ */
