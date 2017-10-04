
#ifndef SRC_NEWMESH_ELEMENTS_PLANE_TRIANGLE3_H_
#define SRC_NEWMESH_ELEMENTS_PLANE_TRIANGLE3_H_

#include "../newelement.h"

namespace espreso {

struct Triangle3 {

	static NewElement create()
	{
		return NewElement(NewElement::TYPE::PLANE, NewElement::CODE::TRIANGLE3, 3, 2, 1);
	}
};

}


#endif /* SRC_NEWMESH_ELEMENTS_PLANE_TRIANGLE3_H_ */
