
#ifndef SRC_NEWMESH_ELEMENTS_VOLUME_PRISMA15_H_
#define SRC_NEWMESH_ELEMENTS_VOLUME_PRISMA15_H_

#include "../element.h"

namespace espreso {

struct Prisma15 {

	static Element create()
	{
		return Element(Element::TYPE::VOLUME, Element::CODE::PRISMA15, 15, 4, 3);
	}
};

}


#endif /* SRC_NEWMESH_ELEMENTS_VOLUME_PRISMA15_H_ */
