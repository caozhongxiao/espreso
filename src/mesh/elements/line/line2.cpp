
#include "line2.h"

#include "../../../basis/containers/serializededata.h"
#include "../../../basis/matrices/denseMatrix.h"

using namespace espreso;

Element Line2::fill(Element e, Element* begin)
{
	size_t GPCount = 2, nodeCount = 2;

	e.N = new std::vector<DenseMatrix>(GPCount, DenseMatrix(1, nodeCount));
	e.dN = new std::vector<DenseMatrix>(GPCount, DenseMatrix(1, nodeCount));
	e.weighFactor = new std::vector<double>(GPCount, 1);

	std::vector<double> s = { 1 / sqrt(3), -1 / sqrt(3) };

	for (unsigned int i = 0; i < GPCount; i++) {
		(*e.N)[i](0, 0) = (1 - s[i]) / 2.0;
		(*e.N)[i](0, 1) = (1 + s[i]) / 2.0;
	}

	for (unsigned int i = 0; i < GPCount; i++) {
		DenseMatrix &m = (*e.dN)[i];

		// dNs - derivation of basis function
		m(0, 0) = -1 / 2.0;
		m(0, 1) =  1 / 2.0;
	}

	return e;
}



