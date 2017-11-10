
#ifndef SRC_MESH_ELEMENTS_POINT_POINT1_H_
#define SRC_MESH_ELEMENTS_POINT_POINT1_H_

#include "../element.h"

namespace espreso {

struct Point1 {

	static Element create()
	{
		return Element(Element::TYPE::POINT, Element::CODE::POINT1, 1, 1, 1);
	}
};

}



#endif /* SRC_MESH_ELEMENTS_POINT_POINT1_H_ */
