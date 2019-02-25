
#include "mesh/elements/element.h"
#include "mesh/mesh.h"

#include "basis/containers/serializededata.h"
#include "basis/matrices/denseMatrix.h"

//using namespace espreso;

namespace espreso {

template<>
void Element::set<Element::CODE::SQUARE8>()
{
	type = Element::TYPE::PLANE;
	code = Element::CODE::SQUARE8;
	nodes = 8;
	coarseNodes = 4;
	nCommonFace = 3;
	nCommonEdge = 2;

	std::vector<Element*> epointers(4, &Mesh::edata[static_cast<int>(Element::CODE::LINE3)]);

	std::vector<int> data = {
		0, 1, 4,
		1, 2, 5,
		2, 3, 6,
		3, 0, 7
	};

	std::vector<int> tringles = {
		0, 4, 7,
		4, 1, 5,
		5, 2, 6,
		6, 3, 7,
		4, 5, 6,
		4, 6, 7
	};

	edges = new serializededata<int, int>(3, data);
	edgepointers = new serializededata<int, Element*>(1, epointers);
	faces = new serializededata<int, int>(3, data);
	facepointers = new serializededata<int, Element*>(1, epointers);
	triangles = new serializededata<int, int>(3, tringles);

	size_t GPCount = 9, nodeCount = 8;

	N = new std::vector<DenseMatrix>(GPCount, DenseMatrix(1, nodeCount));
	dN = new std::vector<DenseMatrix>(GPCount, DenseMatrix(2, nodeCount));
	weighFactor = new std::vector<double>({
			25 / 81.0, 25 / 81.0, 25 / 81.0, 25 / 81.0,
			40 / 81.0, 40 / 81.0, 40 / 81.0, 40 / 81.0,
			64 / 81.0 });

	std::vector< std::vector<double> > st(2, std::vector<double>(GPCount));
	double v = sqrt(0.6);
	st[0] = { -v,  v,  v, -v,  0,  v,  0, -v, 0 };
	st[1] = { -v, -v,  v,  v, -v,  0,  v,  0, 0 };

	for (size_t i = 0; i < GPCount; i++) {
		const std::vector<double> &s = st[0];
		const std::vector<double> &t = st[1];

		(*N)[i](0, 0) = -.25 * (s[i] - 1) * (t[i] - 1) * (s[i] + t[i] + 1);
		(*N)[i](0, 1) =  .25 * (t[i] - 1) * (-s[i] * s[i] + t[i] * s[i] + t[i] + 1);
		(*N)[i](0, 2) =  .25 * (s[i] + 1) * (t[i] + 1) * (s[i] + t[i] - 1);
		(*N)[i](0, 3) =  .25 * (s[i] - 1) * (s[i] - t[i] + 1) * (t[i] + 1);
		(*N)[i](0, 4) =  .5  * (s[i] * s[i] - 1) * (t[i] - 1);
		(*N)[i](0, 5) = -.5  * (s[i] + 1) * (t[i] * t[i] - 1);
		(*N)[i](0, 6) = -.5  * (s[i] * s[i] - 1) * (t[i] + 1);
		(*N)[i](0, 7) =  .5  * (s[i] - 1) * (t[i] * t[i] - 1);
	}

	for (size_t i = 0; i < GPCount; i++) {
		DenseMatrix &m = (*dN)[i];
		const std::vector<double> &s = st[0];
		const std::vector<double> &t = st[1];

		// dNs - derivation of basis function
		m(0, 0) = -((2 * s[i] + t[i]) * (t[i] - 1)) / 4;
		m(0, 1) = -((2 * s[i] - t[i]) * (t[i] - 1)) / 4;
		m(0, 2) = ((2 * s[i] + t[i]) * (t[i] + 1)) / 4;
		m(0, 3) = ((2 * s[i] - t[i]) * (t[i] + 1)) / 4;
		m(0, 4) = s[i] * (t[i] - 1);
		m(0, 5) = 1. / 2 - t[i] * t[i] / 2;
		m(0, 6) = -s[i] * (t[i] + 1);
		m(0, 7) = t[i] * t[i] / 2 - 1. / 2;

		// dNt - derivation of basis function
		m(1, 0) = -((s[i] + 2 * t[i]) * (s[i] - 1)) / 4;
		m(1, 1) = -((s[i] - 2 * t[i]) * (s[i] + 1)) / 4;
		m(1, 2) = ((s[i] + 2 * t[i]) * (s[i] + 1)) / 4;
		m(1, 3) = ((s[i] - 2 * t[i]) * (s[i] - 1)) / 4;
		m(1, 4) = s[i] * s[i] / 2 - 1. / 2;
		m(1, 5) = -t[i] * (s[i] + 1);
		m(1, 6) = 1. / 2 - s[i] * s[i] / 2;
		m(1, 7) = t[i] * (s[i] - 1);
	}
}
}



