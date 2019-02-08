
#include "mesh/elements/element.h"
#include "mesh/mesh.h"

#include "basis/containers/serializededata.h"
#include "basis/matrices/denseMatrix.h"

using namespace espreso;

template<>
void Element::set<Element::CODE::TRIANGLE3>()
{
	type = Element::TYPE::PLANE;
	code = Element::CODE::TRIANGLE3;
	nodes = 3;
	coarseNodes = 3;
	nCommonFace = 2;
	nCommonEdge = 1;

	std::vector<Element*> epointers(3, &Mesh::edata[static_cast<int>(Element::CODE::LINE2)]);

	std::vector<int> data = {
		0, 1,
		1, 2,
		2, 0
	};

	std::vector<int> tringles = {
		0, 1, 2
	};

	edges = new serializededata<int, int>(2, data);
	edgepointers = new serializededata<int, Element*>(1, epointers);
	faces = new serializededata<int, int>(2, data);
	facepointers = new serializededata<int, Element*>(1, epointers);
	triangles = new serializededata<int, int>(3, tringles);

	size_t GPCount = 1, nodeCount = 3;

	N = new std::vector<DenseMatrix>(GPCount, DenseMatrix(1, nodeCount));
	dN = new std::vector<DenseMatrix>(GPCount, DenseMatrix(2, nodeCount));
	weighFactor = new std::vector<double>({ 1 / 2.0 });

	std::vector<double> s = { 1.0 / 3 };
	std::vector<double> t = { 1.0 / 3 };

	for (unsigned int i = 0; i < GPCount; i++) {
		(*N)[i](0, 0) = 1 - s[i] - t[i];
		(*N)[i](0, 1) = s[i];
		(*N)[i](0, 2) = t[i];
	}

	for (unsigned int i = 0; i < GPCount; i++) {
		DenseMatrix &m = (*dN)[i];

		// dNs - derivation of basis function
		m(0, 0) = -1;
		m(0, 1) =  1;
		m(0, 2) =  0;

		// dNt - derivation of basis function
		m(1, 0) = -1;
		m(1, 1) =  0;
		m(1, 2) =  1;
	}
}




