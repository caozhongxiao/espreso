
#ifndef SRC_NEWMESH_ELEMENTS_VOLUME_PRISMA6_H_
#define SRC_NEWMESH_ELEMENTS_VOLUME_PRISMA6_H_

#include "../newelement.h"

namespace espreso {

struct Prisma6 {

	static NewElement create()
	{
		return NewElement(NewElement::TYPE::VOLUME, NewElement::CODE::PRISMA6, 6, 3, 2);
	}
};

}


#endif /* SRC_NEWMESH_ELEMENTS_VOLUME_PRISMA6_H_ */
