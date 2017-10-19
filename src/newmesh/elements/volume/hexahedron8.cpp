
#include "hexahedron8.h"

#include "../../../basis/containers/serializededata.h"

#include "../../../config/ecf/environment.h"

using namespace espreso;

NewElement Hexahedron8::fill(NewElement e, size_t thread)
{
	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<std::vector<int> > data(threads);

	data[thread] = {
		0, 1, 5, 4,
		3, 2, 1, 0,
		4, 5, 6, 7,
		7, 6, 2, 3,
		1, 2, 6, 5,
		3, 0, 4, 7
	};

	e.faces = new serializededata<int, int>(4, data);

	return e;
}




