
#ifndef SRC_MESH_ELEMENTS_VOLUME_PYRAMID5_H_
#define SRC_MESH_ELEMENTS_VOLUME_PYRAMID5_H_

#include "mesh/elements/element.h"

namespace espreso {

struct Pyramid5 {

	static Element fill(Element e, Element* begin);

	static Element create(Element* begin)
	{
		return fill(Element(Element::TYPE::VOLUME, Element::CODE::PYRAMID5, 5, 5, 3, 2), begin);
	}
};

}


#endif /* SRC_MESH_ELEMENTS_VOLUME_PYRAMID5_H_ */
