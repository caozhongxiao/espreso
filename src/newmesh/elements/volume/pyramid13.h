
#ifndef SRC_NEWMESH_ELEMENTS_VOLUME_PYRAMID13_H_
#define SRC_NEWMESH_ELEMENTS_VOLUME_PYRAMID13_H_

#include "../newelement.h"

namespace espreso {

struct Pyramid13 {

	static NewElement create()
	{
		return NewElement(NewElement::TYPE::VOLUME, NewElement::CODE::PYRAMID13, 13, 4, 3);
	}
};

}


#endif /* SRC_NEWMESH_ELEMENTS_VOLUME_PYRAMID13_H_ */
