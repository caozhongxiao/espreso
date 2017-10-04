
#ifndef SRC_NEWMESH_ELEMENTS_VOLUME_PRISMA15_H_
#define SRC_NEWMESH_ELEMENTS_VOLUME_PRISMA15_H_

#include "../newelement.h"

namespace espreso {

struct Prisma15 {

	static NewElement create()
	{
		return NewElement(NewElement::TYPE::VOLUME, NewElement::CODE::PRISMA15, 15, 4, 3);
	}
};

}


#endif /* SRC_NEWMESH_ELEMENTS_VOLUME_PRISMA15_H_ */
