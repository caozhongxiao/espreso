
#include "mesh/elements/element.h"
#include "mesh/mesh.h"

#include "basis/containers/serializededata.h"
#include "basis/matrices/denseMatrix.h"

namespace espreso {

template<>
void Element::set<Element::CODE::PRISMA6>()
{
	type = Element::TYPE::VOLUME;
	code = Element::CODE::PRISMA6;
	nodes = 6;
	coarseNodes = 6;
	nCommonFace = 3;
	nCommonEdge = 2;

	std::vector<Element*> fpointers;
	fpointers.resize(3, &Mesh::edata[static_cast<int>(Element::CODE::SQUARE4)]);
	fpointers.resize(5, &Mesh::edata[static_cast<int>(Element::CODE::TRIANGLE3)]);

	std::vector<int> dist = { 0, 4, 8, 12, 15, 18 };
	std::vector<int> data = {
		0, 1, 4, 3,
		1, 2, 5, 4,
		2, 0, 3, 5,
		2, 1, 0,
		3, 4, 5,
	};

	faces = new serializededata<int, int>(dist, data);
	facepointers = new serializededata<int, Element*>(1, fpointers);

	size_t GPCount = 9, nodeCount = 6;

	N = new std::vector<DenseMatrix>(GPCount, DenseMatrix(1, nodeCount));
	dN = new std::vector<DenseMatrix>(GPCount, DenseMatrix(3, nodeCount));
	weighFactor = new std::vector<double>(GPCount, 1);

	std::vector< std::vector<double> > rst(3, std::vector<double>(GPCount));

	switch (GPCount) {
	case 9: {
		double v1 = 1.0 / 6.0;
		double v2 = 4.0 / 6.0;
		double v3 = sqrt(3.0 / 5.0);
		double v4 = 0.0;
		rst[0] = {  v1,  v2,  v1,  v1,  v2,  v1,  v1,  v2,  v1 };
		rst[1] = {  v1,  v1,  v2,  v1,  v1,  v2,  v1,  v1,  v2 };
		rst[2] = { -v3, -v3, -v3,  v4,  v4,  v4,  v3,  v3,  v3 };
		break;
	}
	default:
		exit(1);
	}

	for (unsigned int i = 0; i < GPCount; i++) {
		DenseMatrix &m = (*N)[i];

		double r = rst[0][i];
		double s = rst[1][i];
		double t = rst[2][i];

		// basis function
		m(0, 0) = 0.5 * ((1.0 - t) * (1.0 - r - s));
		m(0, 1) = 0.5 * ((1.0 - t) * r);
		m(0, 2) = 0.5 * ((1.0 - t) * s);
		m(0, 3) = 0.5 * ((1.0 + t) * (1.0 - r - s));
		m(0, 4) = 0.5 * ((1.0 + t) * r);
		m(0, 5) = 0.5 * ((1.0 + t) * s);
	}

	for (unsigned int i = 0; i < GPCount; i++) {
		///dN contains [dNr, dNs, dNt]
		DenseMatrix &m = (*dN)[i];

		double r = rst[0][i];
		double s = rst[1][i];
		double t = rst[2][i];

		// dNr - derivation of basis function
		m(0, 0) =  t / 2.0 - 1.0 / 2.0;
		m(0, 1) = -t / 2.0 + 1.0 / 2.0;
		m(0, 2) =  0.0;
		m(0, 3) = -t / 2.0 - 1.0 / 2.0;
		m(0, 4) =  t / 2.0 + 1.0 / 2.0;
		m(0, 5) =  0;

		// dNs - derivation of basis function
		m(1, 0) =  t / 2.0 - 1.0 / 2.0;
		m(1, 1) =  0.0;
		m(1, 2) = -t / 2.0 + 1.0 / 2.0;
		m(1, 3) = -t / 2.0 - 1.0 / 2.0;
		m(1, 4) =  0.0;
		m(1, 5) =  t / 2.0 + 1.0 / 2.0;

		// dNt - derivation of basis function
		m(2, 0) =  r / 2.0 + s / 2.0 - 1.0 / 2.0;
		m(2, 1) = -r / 2.0;
		m(2, 2) =          - s / 2.0;
		m(2, 3) = -r / 2.0 - s / 2.0 + 1.0 / 2.0;
		m(2, 4) =  r / 2.0;
		m(2, 5) =            s / 2.0;
	}

	switch (GPCount) {
	case 9: {
		double v1 = 5.0 / 54.0;
		double v2 = 8.0 / 54.0;
		(*weighFactor) = { v1, v1, v1, v2, v2, v2, v1, v1, v1 };
		break;
	}
	default:
		exit(1);
	}
}
}






