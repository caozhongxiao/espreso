
#ifndef SRC_MESH_ELEMENTS_VOLUME_PRISMA6_H_
#define SRC_MESH_ELEMENTS_VOLUME_PRISMA6_H_

#include "../element.h"

namespace espreso {

struct Prisma6 {

	static Element fill(Element e, Element* begin);

	static Element create(Element* begin)
	{
		return fill(Element(Element::TYPE::VOLUME, Element::CODE::PRISMA6, 6, 3, 2), begin);
	}
};

}


#endif /* SRC_MESH_ELEMENTS_VOLUME_PRISMA6_H_ */
