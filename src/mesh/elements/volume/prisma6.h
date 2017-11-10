
#ifndef SRC_MESH_ELEMENTS_VOLUME_PRISMA6_H_
#define SRC_MESH_ELEMENTS_VOLUME_PRISMA6_H_

#include "../element.h"

namespace espreso {

struct Prisma6 {

	static Element create()
	{
		return Element(Element::TYPE::VOLUME, Element::CODE::PRISMA6, 6, 3, 2);
	}
};

}


#endif /* SRC_MESH_ELEMENTS_VOLUME_PRISMA6_H_ */
