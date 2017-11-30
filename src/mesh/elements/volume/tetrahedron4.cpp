
#include "tetrahedron4.h"

#include "../../../basis/containers/serializededata.h"
#include "../../../config/ecf/environment.h"

using namespace espreso;

Element Tetrahedron4::fill(Element e, size_t thread, Element* begin)
{
	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<Element*> facepointers(4, begin + static_cast<int>(Element::CODE::TRIANGLE3));

	std::vector<int> data = {
		0, 1, 3,
		1, 2, 3,
		2, 0, 3,
		2, 1, 0
	};

	e.faces = new serializededata<int, int>(3, data);
	e.facepointers = new serializededata<int, Element*>(1, facepointers);

	return e;
}




