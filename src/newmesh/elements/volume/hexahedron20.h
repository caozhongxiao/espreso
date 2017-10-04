
#ifndef SRC_NEWMESH_ELEMENTS_VOLUME_HEXAHEDRON20_H_
#define SRC_NEWMESH_ELEMENTS_VOLUME_HEXAHEDRON20_H_

#include "../newelement.h"

namespace espreso {

struct Hexahedron20 {

	static NewElement create()
	{
		return NewElement(NewElement::TYPE::VOLUME, NewElement::CODE::HEXA20, 20, 4, 3);
	}
};

}


#endif /* SRC_NEWMESH_ELEMENTS_VOLUME_HEXAHEDRON20_H_ */
