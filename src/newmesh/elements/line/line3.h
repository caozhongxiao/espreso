
#ifndef SRC_NEWMESH_ELEMENTS_LINE_LINE3_H_
#define SRC_NEWMESH_ELEMENTS_LINE_LINE3_H_

#include "../newelement.h"

namespace espreso {

struct Line3 {

	static NewElement create()
	{
		return NewElement(NewElement::TYPE::LINE, NewElement::CODE::LINE3, 3, 1, 1);
	}
};

}


#endif /* SRC_NEWMESH_ELEMENTS_LINE_LINE3_H_ */
