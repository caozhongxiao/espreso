
#ifndef SRC_MESH_ELEMENTS_LINE_LINE2_H_
#define SRC_MESH_ELEMENTS_LINE_LINE2_H_

#include "../element.h"

namespace espreso {

struct Line2 {

	static Element create()
	{
		return Element(Element::TYPE::LINE, Element::CODE::LINE2, 2, 1, 1);
	}
};

}


#endif /* SRC_MESH_ELEMENTS_LINE_LINE2_H_ */
