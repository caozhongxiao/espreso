
#include "triangle3.h"

#include "../../../basis/containers/serializededata.h"
#include "../../../basis/matrices/denseMatrix.h"

using namespace espreso;

Element Triangle3::fill(Element e, Element* begin)
{
	std::vector<Element*> edgepointers(3, begin + static_cast<int>(Element::CODE::LINE2));

	std::vector<int> data = {
		0, 1,
		1, 2,
		2, 0
	};

	e.edges = new serializededata<int, int>(2, data);
	e.edgepointers = new serializededata<int, Element*>(1, edgepointers);
	e.faces = new serializededata<int, int>(2, data);
	e.facepointers = new serializededata<int, Element*>(1, edgepointers);

	size_t GPCount = 1, nodeCount = 3;

	e.N = new std::vector<DenseMatrix>(GPCount, DenseMatrix(1, nodeCount));
	e.dN = new std::vector<DenseMatrix>(GPCount, DenseMatrix(2, nodeCount));
	e.weighFactor = new std::vector<double>({ 1 / 2.0 });

	std::vector<double> s = { 1.0 / 3 };
	std::vector<double> t = { 1.0 / 3 };

	for (unsigned int i = 0; i < GPCount; i++) {
		(*e.N)[i](0, 0) = 1 - s[i] - t[i];
		(*e.N)[i](0, 1) = s[i];
		(*e.N)[i](0, 2) = t[i];
	}

	for (unsigned int i = 0; i < GPCount; i++) {
		DenseMatrix &m = (*e.dN)[i];

		// dNs - derivation of basis function
		m(0, 0) = -1;
		m(0, 1) =  1;
		m(0, 2) =  0;

		// dNt - derivation of basis function
		m(1, 0) = -1;
		m(1, 1) =  0;
		m(1, 2) =  1;
	}

	return e;
}




