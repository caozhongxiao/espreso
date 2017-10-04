
#ifndef SRC_NEWMESH_ELEMENTS_PLANE_SQUARE8_H_
#define SRC_NEWMESH_ELEMENTS_PLANE_SQUARE8_H_

#include "../newelement.h"

namespace espreso {

struct Square8 {

	static NewElement create()
	{
		return NewElement(NewElement::TYPE::PLANE, NewElement::CODE::SQUARE8, 8, 3, 2);
	}
};

}


#endif /* SRC_NEWMESH_ELEMENTS_PLANE_SQUARE8_H_ */
