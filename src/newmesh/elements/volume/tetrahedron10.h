
#ifndef SRC_NEWMESH_ELEMENTS_VOLUME_TETRAHEDRON10_H_
#define SRC_NEWMESH_ELEMENTS_VOLUME_TETRAHEDRON10_H_

#include "../newelement.h"

namespace espreso {

struct Tetrahedron10 {

	static NewElement create()
	{
		return NewElement(NewElement::TYPE::VOLUME, NewElement::CODE::TETRA10, 10, 4, 3);
	}
};

}

#endif /* SRC_NEWMESH_ELEMENTS_VOLUME_TETRAHEDRON10_H_ */
