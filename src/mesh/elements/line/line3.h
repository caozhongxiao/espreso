
#ifndef SRC_MESH_ELEMENTS_LINE_LINE3_H_
#define SRC_MESH_ELEMENTS_LINE_LINE3_H_

#include "mesh/elements/element.h"

namespace espreso {

struct Line3 {

	static Element fill(Element e, Element* begin);

	static Element create(Element* begin)
	{
		return fill(Element(Element::TYPE::LINE, Element::CODE::LINE3, 3, 2, 1, 1), begin);
	}
};

}


#endif /* SRC_MESH_ELEMENTS_LINE_LINE3_H_ */
