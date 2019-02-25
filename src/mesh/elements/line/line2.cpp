
#include "mesh/elements/element.h"

#include "basis/containers/serializededata.h"
#include "basis/matrices/denseMatrix.h"

//using namespace espreso;

namespace espreso {
template<>
void Element::set<Element::CODE::LINE2>()
{
	type = Element::TYPE::LINE;
	code = Element::CODE::LINE2;
	nodes = 2;
	coarseNodes = 2;
	nCommonFace = 1;
	nCommonEdge = 1;

	size_t GPCount = 2, nodeCount = 2;

	N = new std::vector<DenseMatrix>(GPCount, DenseMatrix(1, nodeCount));
	dN = new std::vector<DenseMatrix>(GPCount, DenseMatrix(1, nodeCount));
	weighFactor = new std::vector<double>(GPCount, 1);

	std::vector<double> s = { 1 / sqrt(3), -1 / sqrt(3) };

	for (unsigned int i = 0; i < GPCount; i++) {
		(*N)[i](0, 0) = (1 - s[i]) / 2.0;
		(*N)[i](0, 1) = (1 + s[i]) / 2.0;
	}

	for (unsigned int i = 0; i < GPCount; i++) {
		DenseMatrix &m = (*dN)[i];

		// dNs - derivation of basis function
		m(0, 0) = -1 / 2.0;
		m(0, 1) =  1 / 2.0;
	}
}
}


