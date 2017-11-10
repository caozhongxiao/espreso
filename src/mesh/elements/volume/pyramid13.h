
#ifndef SRC_MESH_ELEMENTS_VOLUME_PYRAMID13_H_
#define SRC_MESH_ELEMENTS_VOLUME_PYRAMID13_H_

#include "../element.h"

namespace espreso {

struct Pyramid13 {

	static Element create()
	{
		return Element(Element::TYPE::VOLUME, Element::CODE::PYRAMID13, 13, 4, 3);
	}
};

}


#endif /* SRC_MESH_ELEMENTS_VOLUME_PYRAMID13_H_ */
