
#ifndef SRC_MESH_ELEMENTS_LINE_LINE3_H_
#define SRC_MESH_ELEMENTS_LINE_LINE3_H_

#include "../element.h"

namespace espreso {

struct Line3 {

	static Element create()
	{
		return Element(Element::TYPE::LINE, Element::CODE::LINE3, 3, 1, 1);
	}
};

}


#endif /* SRC_MESH_ELEMENTS_LINE_LINE3_H_ */
