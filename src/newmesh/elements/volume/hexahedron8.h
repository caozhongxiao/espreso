
#ifndef SRC_NEWMESH_ELEMENTS_VOLUME_HEXAHEDRON8_H_
#define SRC_NEWMESH_ELEMENTS_VOLUME_HEXAHEDRON8_H_

#include "../newelement.h"

namespace espreso {

struct Hexahedron8 {

	static NewElement fill(NewElement e, size_t thread);

	static NewElement create(size_t thread)
	{
		return fill(NewElement(NewElement::TYPE::VOLUME, NewElement::CODE::HEXA8, 8, 3, 2), thread);
	}

};

}


#endif /* SRC_NEWMESH_ELEMENTS_VOLUME_HEXAHEDRON8_H_ */
