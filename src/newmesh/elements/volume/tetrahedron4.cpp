
#include "tetrahedron4.h"

#include "../../../basis/containers/serializededata.h"
#include "../../../config/ecf/environment.h"

using namespace espreso;

NewElement Tetrahedron4::fill(NewElement e, size_t thread, NewElement* begin)
{
	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<NewElement*> facepointers(4, begin + static_cast<int>(NewElement::CODE::TRIANGLE3));

	std::vector<int> data = {
		0, 1, 3,
		1, 2, 3,
		2, 0, 3,
		2, 1, 0
	};

	e.faces = new serializededata<int, int>(3, { thread, threads, data });
	e.facepointers = new serializededata<int, NewElement*>(1, { thread, threads, facepointers });

	return e;
}




