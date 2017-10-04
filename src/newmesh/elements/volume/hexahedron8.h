
#ifndef SRC_NEWMESH_ELEMENTS_VOLUME_HEXAHEDRON8_H_
#define SRC_NEWMESH_ELEMENTS_VOLUME_HEXAHEDRON8_H_

#include "../newelement.h"

namespace espreso {

struct Hexahedron8 {

	static NewElement create()
	{
		return NewElement(NewElement::TYPE::VOLUME, NewElement::CODE::HEXA8, 8, 3, 2);
	}
};

}


#endif /* SRC_NEWMESH_ELEMENTS_VOLUME_HEXAHEDRON8_H_ */
