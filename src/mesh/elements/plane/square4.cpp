
#include "square4.h"

#include "../../../basis/containers/serializededata.h"
#include "../../../basis/matrices/denseMatrix.h"

using namespace espreso;

Element Square4::fill(Element e, Element* begin)
{
	std::vector<Element*> edgepointers(4, begin + static_cast<int>(Element::CODE::LINE2));

	std::vector<int> data = {
		0, 1,
		1, 2,
		2, 3,
		3, 0
	};

	e.edges = new serializededata<int, int>(2, data);
	e.edgepointers = new serializededata<int, Element*>(1, edgepointers);
	e.faces = new serializededata<int, int>(2, data);
	e.facepointers = new serializededata<int, Element*>(1, edgepointers);

	size_t GPCount = 4, nodeCount = 4;

	e.N = new std::vector<DenseMatrix>(GPCount, DenseMatrix(1, nodeCount));
	e.dN = new std::vector<DenseMatrix>(GPCount, DenseMatrix(2, nodeCount));
	e.weighFactor = new std::vector<double>(GPCount, 1);

	double CsQ_scale = 0.577350269189626;

	std::vector<double> s = { -CsQ_scale,  CsQ_scale,  CsQ_scale, -CsQ_scale };
	std::vector<double> t = { -CsQ_scale, -CsQ_scale,  CsQ_scale,  CsQ_scale };

	for (unsigned int i = 0; i < GPCount; i++) {
		(*e.N)[i](0, 0) = 0.25 * (1 - s[i]) * (1 - t[i]);
		(*e.N)[i](0, 1) = 0.25 * (s[i] + 1) * (1 - t[i]);
		(*e.N)[i](0, 2) = 0.25 * (s[i] + 1) * (t[i] + 1);
		(*e.N)[i](0, 3) = 0.25 * (1 - s[i]) * (t[i] + 1);
	}

	for (unsigned int i = 0; i < GPCount; i++) {
		// dNs - derivation of basis function
		(*e.dN)[i](0, 0) = 0.25 * ( t[i] - 1);
		(*e.dN)[i](0, 1) = 0.25 * (-t[i] + 1);
		(*e.dN)[i](0, 2) = 0.25 * ( t[i] + 1);
		(*e.dN)[i](0, 3) = 0.25 * (-t[i] - 1);

		// dNt - derivation of basis function
		(*e.dN)[i](1, 0) = 0.25 * ( s[i] - 1);
		(*e.dN)[i](1, 1) = 0.25 * (-s[i] - 1);
		(*e.dN)[i](1, 2) = 0.25 * ( s[i] + 1);
		(*e.dN)[i](1, 3) = 0.25 * (-s[i] + 1);
	}

	return e;
}



