
#ifndef SRC_NEWMESH_ELEMENTS_LINE_LINE2_H_
#define SRC_NEWMESH_ELEMENTS_LINE_LINE2_H_

#include "../newelement.h"

namespace espreso {

struct Line2 {

	static NewElement create()
	{
		return NewElement(NewElement::TYPE::LINE, NewElement::CODE::LINE2, 2, 1, 1);
	}
};

}


#endif /* SRC_NEWMESH_ELEMENTS_LINE_LINE2_H_ */
