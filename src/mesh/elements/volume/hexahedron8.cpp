
#include "hexahedron8.h"

#include "../../../basis/containers/serializededata.h"
#include "../../../config/ecf/environment.h"

using namespace espreso;

Element Hexahedron8::fill(Element e, size_t thread, Element* begin)
{
	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<Element*> facepointers(6, begin + static_cast<int>(Element::CODE::SQUARE4));

	std::vector<int> data = {
		0, 1, 5, 4,
		3, 2, 1, 0,
		4, 5, 6, 7,
		7, 6, 2, 3,
		1, 2, 6, 5,
		3, 0, 4, 7
	};

	e.faces = new serializededata<int, int>(4, { thread, threads, data });
	e.facepointers = new serializededata<int, Element*>(1, { thread, threads, facepointers });

	return e;
}





