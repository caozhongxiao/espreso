
#ifndef SRC_NEWMESH_ELEMENTS_PLANE_SQUARE4_H_
#define SRC_NEWMESH_ELEMENTS_PLANE_SQUARE4_H_

#include "../newelement.h"

namespace espreso {

struct Square4 {

	static NewElement create()
	{
		return NewElement(NewElement::TYPE::PLANE, NewElement::CODE::SQUARE4, 4, 2, 1);
	}
};

}


#endif /* SRC_NEWMESH_ELEMENTS_PLANE_SQUARE4_H_ */
