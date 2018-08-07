
#ifndef SRC_MESH_ELEMENTS_LINE_LINE2_H_
#define SRC_MESH_ELEMENTS_LINE_LINE2_H_

#include "../element.h"

namespace espreso {

struct Line2 {

	static Element fill(Element e, Element* begin);

	static Element create(Element* begin)
	{
		return fill(Element(Element::TYPE::LINE, Element::CODE::LINE2, 2, 2, 1, 1), begin);
	}
};

}


#endif /* SRC_MESH_ELEMENTS_LINE_LINE2_H_ */
