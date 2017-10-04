
#ifndef SRC_NEWMESH_ELEMENTS_VOLUME_TETRAHEDRON4_H_
#define SRC_NEWMESH_ELEMENTS_VOLUME_TETRAHEDRON4_H_

#include "../newelement.h"

namespace espreso {

struct Tetrahedron4 {

	static NewElement create()
	{
		return NewElement(NewElement::TYPE::VOLUME, NewElement::CODE::TETRA4, 4, 3, 2);
	}
};

}


#endif /* SRC_NEWMESH_ELEMENTS_VOLUME_TETRAHEDRON4_H_ */
