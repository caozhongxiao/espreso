
#ifndef SRC_NEWMESH_ELEMENTS_VOLUME_PYRAMID5_H_
#define SRC_NEWMESH_ELEMENTS_VOLUME_PYRAMID5_H_

#include "../newelement.h"

namespace espreso {

struct Pyramid5 {

	static NewElement create()
	{
		return NewElement(NewElement::TYPE::VOLUME, NewElement::CODE::PYRAMID5, 5, 3, 2);
	}
};

}


#endif /* SRC_NEWMESH_ELEMENTS_VOLUME_PYRAMID5_H_ */
