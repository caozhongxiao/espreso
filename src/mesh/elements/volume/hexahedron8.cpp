
#include "mesh/elements/element.h"
#include "mesh/mesh.h"

#include "basis/containers/serializededata.h"
#include "basis/matrices/denseMatrix.h"

namespace espreso {

template<>
void Element::set<Element::CODE::HEXA8>()
{
	type = Element::TYPE::VOLUME;
	code = Element::CODE::HEXA8;
	nodes = 8;
	coarseNodes = 8;
	nCommonFace = 3;
	nCommonEdge = 2;

	std::vector<Element*> fpointers(6, &Mesh::edata[static_cast<int>(Element::CODE::SQUARE4)]);

	std::vector<int> data = {
		0, 1, 5, 4,
		3, 2, 1, 0,
		4, 5, 6, 7,
		7, 6, 2, 3,
		1, 2, 6, 5,
		3, 0, 4, 7
	};

	faces = new serializededata<int, int>(4, data);
	facepointers = new serializededata<int, Element*>(1, fpointers);

	size_t GPCount = 8, nodeCount = 8;

	N = new std::vector<DenseMatrix>(GPCount, DenseMatrix(1, nodeCount));
	dN = new std::vector<DenseMatrix>(GPCount, DenseMatrix(3, nodeCount));
	weighFactor = new std::vector<double>(GPCount, 1);

	double CsQ_scale = 0.577350269189626;

	for (unsigned int i = 0; i < GPCount; i++) {
		double r = (i & 4) ? CsQ_scale : -CsQ_scale;
		double s = (i & 2) ? CsQ_scale : -CsQ_scale;
		double t = (i & 1) ? CsQ_scale : -CsQ_scale;

		// basis function
		(*N)[i](0, 0) = 0.125 * (1 - r) * (1 - s) * (1 - t);
		(*N)[i](0, 1) = 0.125 * (r + 1) * (1 - s) * (1 - t);
		(*N)[i](0, 2) = 0.125 * (r + 1) * (s + 1) * (1 - t);
		(*N)[i](0, 3) = 0.125 * (1 - r) * (s + 1) * (1 - t);
		(*N)[i](0, 4) = 0.125 * (1 - r) * (1 - s) * (t + 1);
		(*N)[i](0, 5) = 0.125 * (r + 1) * (1 - s) * (t + 1);
		(*N)[i](0, 6) = 0.125 * (r + 1) * (s + 1) * (t + 1);
		(*N)[i](0, 7) = 0.125 * (1 - r) * (s + 1) * (t + 1);
	}

	for (unsigned int i = 0; i < GPCount; i++) {
		double r = (i & 4) ? CsQ_scale : -CsQ_scale;
		double s = (i & 2) ? CsQ_scale : -CsQ_scale;
		double t = (i & 1) ? CsQ_scale : -CsQ_scale;

		// dNr - derivation of basis function
		(*dN)[i](0, 0) = 0.125 * (-(1 - s) * (1 - t));
		(*dN)[i](0, 1) = 0.125 * ( (1 - s) * (1 - t));
		(*dN)[i](0, 2) = 0.125 * ( (1 + s) * (1 - t));
		(*dN)[i](0, 3) = 0.125 * (-(1 + s) * (1 - t));
		(*dN)[i](0, 4) = 0.125 * (-(1 - s) * (1 + t));
		(*dN)[i](0, 5) = 0.125 * ( (1 - s) * (1 + t));
		(*dN)[i](0, 6) = 0.125 * ( (1 + s) * (1 + t));
		(*dN)[i](0, 7) = 0.125 * (-(1 + s) * (1 + t));

		// dNs - derivation of basis function
		(*dN)[i](1, 0)  = 0.125 * (-(1 - r) * (1 - t));
		(*dN)[i](1, 1)  = 0.125 * (-(1 + r) * (1 - t));
		(*dN)[i](1, 2) = 0.125 * ( (1 + r) * (1 - t));
		(*dN)[i](1, 3) = 0.125 * ( (1 - r) * (1 - t));
		(*dN)[i](1, 4) = 0.125 * (-(1 - r) * (1 + t));
		(*dN)[i](1, 5) = 0.125 * (-(1 + r) * (1 + t));
		(*dN)[i](1, 6) = 0.125 * ( (1 + r) * (1 + t));
		(*dN)[i](1, 7) = 0.125 * ( (1 - r) * (1 + t));

		// dNt - derivation of basis function
		(*dN)[i](2, 0) = 0.125 * (-(1 - r) * (1 - s));
		(*dN)[i](2, 1) = 0.125 * (-(1 + r) * (1 - s));
		(*dN)[i](2, 2) = 0.125 * (-(1 + r) * (1 + s));
		(*dN)[i](2, 3) = 0.125 * (-(1 - r) * (1 + s));
		(*dN)[i](2, 4) = 0.125 * ( (1 - r) * (1 - s));
		(*dN)[i](2, 5) = 0.125 * ( (1 + r) * (1 - s));
		(*dN)[i](2, 6) = 0.125 * ( (1 + r) * (1 + s));
		(*dN)[i](2, 7) = 0.125 * ( (1 - r) * (1 + s));
	}
}
}




