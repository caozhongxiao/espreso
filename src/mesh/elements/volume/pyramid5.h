
#ifndef SRC_MESH_ELEMENTS_VOLUME_PYRAMID5_H_
#define SRC_MESH_ELEMENTS_VOLUME_PYRAMID5_H_

#include "../element.h"

namespace espreso {

struct Pyramid5 {

	static Element create()
	{
		return Element(Element::TYPE::VOLUME, Element::CODE::PYRAMID5, 5, 3, 2);
	}
};

}


#endif /* SRC_MESH_ELEMENTS_VOLUME_PYRAMID5_H_ */
