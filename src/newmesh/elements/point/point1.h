
#ifndef SRC_NEWMESH_ELEMENTS_POINT_POINT1_H_
#define SRC_NEWMESH_ELEMENTS_POINT_POINT1_H_

#include "../newelement.h"

namespace espreso {

struct Point1 {

	static NewElement create()
	{
		return NewElement(NewElement::TYPE::POINT, NewElement::CODE::POINT1, 1, 1, 1);
	}
};

}



#endif /* SRC_NEWMESH_ELEMENTS_POINT_POINT1_H_ */
