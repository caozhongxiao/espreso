
#include "mesh/elements/element.h"
#include "mesh/mesh.h"

#include "basis/containers/serializededata.h"
#include "basis/matrices/denseMatrix.h"

using namespace espreso;

template<>
void Element::set<Element::CODE::SQUARE4>()
{
	type = Element::TYPE::PLANE;
	code = Element::CODE::SQUARE4;
	nodes = 4;
	coarseNodes = 4;
	nCommonFace = 3;
	nCommonEdge = 2;

	std::vector<Element*> epointers(4, &Mesh::edata[static_cast<int>(Element::CODE::LINE2)]);

	std::vector<int> lines = {
		0, 1,
		1, 2,
		2, 3,
		3, 0
	};

	std::vector<int> tringles = {
		0, 1, 2,
		0, 2, 3
	};

	edges = new serializededata<int, int>(2, lines);
	edgepointers = new serializededata<int, Element*>(1, epointers);
	faces = new serializededata<int, int>(2, lines);
	facepointers = new serializededata<int, Element*>(1, epointers);
	triangles = new serializededata<int, int>(3, tringles);

	size_t GPCount = 4, nodeCount = 4;

	N = new std::vector<DenseMatrix>(GPCount, DenseMatrix(1, nodeCount));
	dN = new std::vector<DenseMatrix>(GPCount, DenseMatrix(2, nodeCount));
	weighFactor = new std::vector<double>(GPCount, 1);

	double CsQ_scale = 0.577350269189626;

	std::vector<double> s = { -CsQ_scale,  CsQ_scale,  CsQ_scale, -CsQ_scale };
	std::vector<double> t = { -CsQ_scale, -CsQ_scale,  CsQ_scale,  CsQ_scale };

	for (unsigned int i = 0; i < GPCount; i++) {
		(*N)[i](0, 0) = 0.25 * (1 - s[i]) * (1 - t[i]);
		(*N)[i](0, 1) = 0.25 * (s[i] + 1) * (1 - t[i]);
		(*N)[i](0, 2) = 0.25 * (s[i] + 1) * (t[i] + 1);
		(*N)[i](0, 3) = 0.25 * (1 - s[i]) * (t[i] + 1);
	}

	for (unsigned int i = 0; i < GPCount; i++) {
		// dNs - derivation of basis function
		(*dN)[i](0, 0) = 0.25 * ( t[i] - 1);
		(*dN)[i](0, 1) = 0.25 * (-t[i] + 1);
		(*dN)[i](0, 2) = 0.25 * ( t[i] + 1);
		(*dN)[i](0, 3) = 0.25 * (-t[i] - 1);

		// dNt - derivation of basis function
		(*dN)[i](1, 0) = 0.25 * ( s[i] - 1);
		(*dN)[i](1, 1) = 0.25 * (-s[i] - 1);
		(*dN)[i](1, 2) = 0.25 * ( s[i] + 1);
		(*dN)[i](1, 3) = 0.25 * (-s[i] + 1);
	}
}



