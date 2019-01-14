
#ifndef SRC_MESH_ELEMENTS_VOLUME_PRISMA15_H_
#define SRC_MESH_ELEMENTS_VOLUME_PRISMA15_H_

#include "mesh/elements/element.h"

namespace espreso {

struct Prisma15 {

	static Element fill(Element e, Element* begin);

	static Element create(Element* begin)
	{
		return fill(Element(Element::TYPE::VOLUME, Element::CODE::PRISMA15, 15, 6, 4, 3), begin);
	}
};

}


#endif /* SRC_MESH_ELEMENTS_VOLUME_PRISMA15_H_ */
