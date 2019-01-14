
#include "line3.h"

#include "basis/containers/serializededata.h"
#include "basis/matrices/denseMatrix.h"

using namespace espreso;

Element Line3::fill(Element e, Element* begin)
{
	size_t GPCount = 3, nodeCount = 3;

	e.N = new std::vector<DenseMatrix>(GPCount, DenseMatrix(1, nodeCount));
	e.dN = new std::vector<DenseMatrix>(GPCount, DenseMatrix(1, nodeCount));
	e.weighFactor = new std::vector<double>({ 5 / 9.0, 8 / 9.0, 5 / 9.0 });

	std::vector<double> s = { -sqrt(3 / 5.0), 0, sqrt(3 / 5.0) };

	for (size_t i = 0; i < GPCount; i++) {
		(*e.N)[i](0, 0) = (1 / 2.0) * (s[i] - 1) * s[i];
		(*e.N)[i](0, 1) = 1 - s[i] * s[i];
		(*e.N)[i](0, 2) = (1 / 2.0) * (s[i] + 1) * s[i];
	}

	for (unsigned int i = 0; i < GPCount; i++) {
		DenseMatrix &m = (*e.dN)[i];

		// dNs - derivation of basis function
		m(0, 0) = s[i] - 1 / 2.0;
		m(0, 1) = -2 * s[i];
		m(0, 2) = s[i] + 1 / 2.0;
	}

	return e;
}



