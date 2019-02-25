
#include "mesh/elements/element.h"

#include "basis/containers/serializededata.h"
#include "basis/matrices/denseMatrix.h"

namespace espreso {

template<>
void Element::set<Element::CODE::LINE3>()
{
	type = Element::TYPE::LINE;
	code = Element::CODE::LINE3;
	nodes = 3;
	coarseNodes = 2;
	nCommonFace = 1;
	nCommonEdge = 1;

	size_t GPCount = 3, nodeCount = 3;

	N = new std::vector<DenseMatrix>(GPCount, DenseMatrix(1, nodeCount));
	dN = new std::vector<DenseMatrix>(GPCount, DenseMatrix(1, nodeCount));
	weighFactor = new std::vector<double>({ 5 / 9.0, 8 / 9.0, 5 / 9.0 });

	std::vector<double> s = { -sqrt(3 / 5.0), 0, sqrt(3 / 5.0) };

	for (size_t i = 0; i < GPCount; i++) {
		(*N)[i](0, 0) = (1 / 2.0) * (s[i] - 1) * s[i];
		(*N)[i](0, 1) = 1 - s[i] * s[i];
		(*N)[i](0, 2) = (1 / 2.0) * (s[i] + 1) * s[i];
	}

	for (unsigned int i = 0; i < GPCount; i++) {
		DenseMatrix &m = (*dN)[i];

		// dNs - derivation of basis function
		m(0, 0) = s[i] - 1 / 2.0;
		m(0, 1) = -2 * s[i];
		m(0, 2) = s[i] + 1 / 2.0;
	}
}
}


