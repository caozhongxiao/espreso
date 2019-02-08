
#include "mesh/elements/element.h"
#include "mesh/mesh.h"

#include "basis/containers/serializededata.h"
#include "basis/matrices/denseMatrix.h"

using namespace espreso;

template<>
void Element::set<Element::CODE::TRIANGLE6>()
{
	type = Element::TYPE::PLANE;
	code = Element::CODE::TRIANGLE6;
	nodes = 6;
	coarseNodes = 3;
	nCommonFace = 3;
	nCommonEdge = 2;

	std::vector<Element*> epointers(3, &Mesh::edata[static_cast<int>(Element::CODE::LINE3)]);

	std::vector<int> data = {
		0, 1, 3,
		1, 2, 4,
		2, 0, 5
	};

	std::vector<int> tringles = {
		0, 3, 5,
		3, 1, 4,
		4, 2, 5,
		5, 0, 3,
		3, 4, 5
	};

	edges = new serializededata<int, int>(3, data);
	edgepointers = new serializededata<int, Element*>(1, epointers);
	faces = new serializededata<int, int>(3, data);
	facepointers = new serializededata<int, Element*>(1, epointers);
	triangles = new serializededata<int, int>(3, tringles);

	size_t GPCount = 6, nodeCount = 6;

	N = new std::vector<DenseMatrix>(GPCount, DenseMatrix(1, nodeCount));
	dN = new std::vector<DenseMatrix>(GPCount, DenseMatrix(2, nodeCount));
	weighFactor = new std::vector<double>(
		{ 0.109951743655322 / 2.0, 0.109951743655322 / 2.0, 0.109951743655322 / 2.0,
		  0.223381589678011 / 2.0, 0.223381589678011 / 2.0, 0.223381589678011 / 2.0 });

	const std::vector<double> s = {
		0.091576213509771,
		0.816847572980459,
		0.091576213509771,
		0.445948490915965,
		0.108103018168070,
		0.445948490915965
	};
	const std::vector<double> t = {
		0.091576213509771,
		0.091576213509771,
		0.816847572980459,
		0.445948490915965,
		0.445948490915965,
		0.108103018168070
	};

	for (size_t i = 0; i < GPCount; i++) {
		(*N)[i](0, 0) = (1.0 - s[i] - t[i]) * (1.0 - 2.0 * (s[i] + t[i]));
		(*N)[i](0, 1) = -(s[i]) * (1.0 - 2.0 * s[i]);
		(*N)[i](0, 2) = -(t[i]) * (1.0 - 2.0 * t[i]);
		(*N)[i](0, 3) = 4.0 * (s[i]) * (1.0 - s[i] - t[i]);
		(*N)[i](0, 4) = 4.0 * (s[i]) * (t[i]);
		(*N)[i](0, 5) = 4.0 * (t[i]) * (1.0 - s[i] - t[i]);
	}

	for (size_t i = 0; i < GPCount; i++) {
		///dN contains [dNs, dNt]
		DenseMatrix &m = (*dN)[i];

		// dNs - derivation of basis function
		m(0, 0) = -3.0 + 4.0 * s[i] + 4.0 * t[i];
		m(0, 1) = -1.0 + 4.0 * s[i];
		m(0, 2) = 0.0;
		m(0, 3) = 4.0 - 8.0 * s[i] - 4.0 * t[i];
		m(0, 4) = 4.0 * t[i];
		m(0, 5) = -4.0 * t[i];

		// dNt - derivation of basis function
		m(1, 0) = -3.0 + 4.0 * s[i] + 4.0 * t[i];
		m(1, 1) = 0.0;
		m(1, 2) = -1.0 + 4.0 * t[i];
		m(1, 3) = -4.0 * s[i];
		m(1, 4) = 4.0 * s[i];
		m(1, 5) = 4.0 - 4.0 * s[i] - 8.0 * t[i];
	}
}
