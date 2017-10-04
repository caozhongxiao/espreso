
#ifndef SRC_NEWMESH_ELEMENTS_PLANE_TRIANGLE6_H_
#define SRC_NEWMESH_ELEMENTS_PLANE_TRIANGLE6_H_

#include "../newelement.h"

namespace espreso {

struct Triangle6 {

	static NewElement create()
	{
		return NewElement(NewElement::TYPE::PLANE, NewElement::CODE::TRIANGLE6, 6, 3, 2);
	}
};

}


#endif /* SRC_NEWMESH_ELEMENTS_PLANE_TRIANGLE6_H_ */
