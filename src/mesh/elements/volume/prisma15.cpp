
#include "mesh/elements/element.h"
#include "mesh/mesh.h"

#include "basis/containers/serializededata.h"
#include "basis/matrices/denseMatrix.h"

//using namespace espreso;

namespace espreso {
template<>
void Element::set<Element::CODE::PRISMA15>()
{
	type = Element::TYPE::VOLUME;
	code = Element::CODE::PRISMA15;
	nodes = 15;
	coarseNodes = 6;
	nCommonFace = 4;


	std::vector<Element*> fpointers;
	fpointers.resize(3, &Mesh::edata[static_cast<int>(Element::CODE::SQUARE8)]);
	fpointers.resize(5, &Mesh::edata[static_cast<int>(Element::CODE::TRIANGLE6)]);

	std::vector<int> dist = { 0, 8, 16, 24, 30, 36 };
	std::vector<int> data = {
		0, 1, 4, 3,  6, 13,  9, 12,
		1, 2, 5, 4,  7, 14, 10, 13,
		2, 0, 3, 5,  8, 12, 11, 14,
		2, 1, 0, 7,  6,  8,
		3, 4, 5, 9, 10, 11
	};

	faces = new serializededata<int, int>(dist, data);
	facepointers = new serializededata<int, Element*>(1, fpointers);

	size_t GPCount = 9, nodeCount = 15;

	N = new std::vector<DenseMatrix>(GPCount, DenseMatrix(1, nodeCount));
	dN = new std::vector<DenseMatrix>(GPCount, DenseMatrix(3, nodeCount));
	weighFactor = new std::vector<double>();

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
		m(0, 0) = -(1.0 - r - s) * (1.0 - t) * (2.0 * r + 2.0 * s + t) / 2.0;
		m(0, 1) = r * (1.0 - t) * (2.0 * r - t - 2.0) / 2.0;
		m(0, 2) = s * (1.0 - t) * (2.0 * s - t - 2.0) / 2.0;
		m(0, 3) = -(1.0 - r - s) * (1.0 + t) * (2.0 * r + 2.0 * s - t) / 2.0;
		m(0, 4) = r * (t + 1.0) * (2.0 * r + t - 2.0) / 2.0;
		m(0, 5) = s * (t + 1.0) * (2.0 * s + t - 2.0) / 2.0;
		m(0, 6) = 2.0 * r * (1.0 - r - s) * (1.0 - t);
		m(0, 7) = 2.0 * r * s * (1.0 - t);
		m(0, 8) = 2.0 * s * (1.0 - r - s) * (1.0 - t);
		m(0, 9) = 2.0 * r * (1.0 - r - s) * (1.0 + t);
		m(0, 10) = 2.0 * r * s * (1.0 + t);
		m(0, 11) = 2.0 * s * (1.0 - r - s) * (1.0 + t);
		m(0, 12) = (1.0 - r - s) * (1.0 - pow(t, 2.0));
		m(0, 13) = r * (1.0 - pow(t, 2.0));
		m(0, 14) = s * (1.0 - pow(t, 2.0));
	}

	for (unsigned int i = 0; i < GPCount; i++) {
		///dN contains [dNr, dNs, dNt]
		DenseMatrix &m = (*dN)[i];

		double r = rst[0][i];
		double s = rst[1][i];
		double t = rst[2][i];

		// dNr - derivation of basis function
		m(0, 0) = -(t - 1.0) * (r + s - 1.0) - ((t - 1.0) * (2.0 * r + 2.0 * s + t)) / 2.0;
		m(0, 1) = ((t - 1.0) * (t - 2.0 * r + 2.0)) / 2.0 - r * (t - 1.0);
		m(0, 2) = 0.0;
		m(0, 3) = ((t + 1.0) * (2.0 * r + 2.0 * s - t)) / 2.0 + (t + 1.0) * (r + s - 1.0);
		m(0, 4) = r * (t + 1.0) + ((t + 1.0) * (2.0 * r + t - 2.0)) / 2.0;
		m(0, 5) = 0.0;
		m(0, 6) = 2.0 * (t - 1.0) * (r + s - 1.0) + 2.0 * r * (t - 1.0);
		m(0, 7) = (-2.0) * s * (t - 1.0);
		m(0, 8) = 2.0 * s * (t - 1.0);
		m(0, 9) = -2.0 * (t + 1.0) * (r + s - 1.0) - 2.0 * r * (t + 1.0);
		m(0, 10) = 2.0 * s * (t + 1.0);
		m(0, 11) = -2.0 * s * (t + 1.0);
		m(0, 12) = pow(t, 2.0) - 1.0;
		m(0, 13) = 1.0 - pow(t, 2.0);
		m(0, 14) = 0.0;

		// dNs - derivation of basis function
		m(1, 0) = -(t - 1.0) * (r + s - 1.0) - ((t - 1.0) * (2.0 * r + 2.0 * s + t)) / 2.0;
		m(1, 1) = 0.0;
		m(1, 2) = ((t - 1.0) * (t - 2.0 * s + 2.0)) / 2.0 - s * (t - 1.0);
		m(1, 3) = ((t + 1.0) * (2.0 * r + 2.0 * s - t)) / 2.0 + (t + 1.0) * (r + s - 1.0);
		m(1, 4) = 0.0;
		m(1, 5) = s * (t + 1.0) + ((t + 1.0) * (2.0 * s + t - 2.0)) / 2.0;
		m(1, 6) = 2.0 * r * (t - 1.0);
		m(1, 7) = (-2.0) * r * (t - 1.0);
		m(1, 8) = 2.0 * (t - 1.0) * (r + s - 1.0) + 2.0 * s * (t - 1.0);
		m(1, 9) = (-2.0) * r * (t + 1.0);
		m(1, 10) = 2.0 * r * (t + 1.0);
		m(1, 11) = -2.0 * (t + 1.0) * (r + s - 1.0) - 2.0 * s * (t + 1.0);
		m(1, 12) = pow(t, 2.0) - 1.0;
		m(1, 13) = 0.0;
		m(1, 14) = 1.0 - pow(t, 2.0);

		// dNt - derivation of basis function
		m(2, 0) = -((r + s - 1.0) * (2.0 * r + 2.0 * s + t)) / 2.0 - ((t - 1.0) * (r + s - 1.0)) / 2.0;
		m(2, 1) = (r * (t - 2.0 * r + 2.0)) / 2.0 + (r * (t - 1.0)) / 2.0;
		m(2, 2) = (s * (t - 2.0 * s + 2.0)) / 2.0 + (s * (t - 1.0)) / 2.0;
		m(2, 3) = ((r + s - 1.0) * (2.0 * r + 2.0 * s - t)) / 2.0 - ((t + 1.0) * (r + s - 1.0)) / 2.0;
		m(2, 4) = (r * (2.0 * r + t - 2.0)) / 2.0 + (r * (t + 1.0)) / 2.0;
		m(2, 5) = (s * (2.0 * s + t - 2.0)) / 2.0 + (s * (t + 1.0)) / 2.0;
		m(2, 6) = 2.0 * r * (r + s - 1.0);
		m(2, 7) = (-2.0) * r * s;
		m(2, 8) = 2.0 * s * (r + s - 1.0);
		m(2, 9) = (-2.0) * r * (r + s - 1.0);
		m(2, 10) = 2.0 * r * s;
		m(2, 11) = (-2.0) * s * (r + s - 1.0);
		m(2, 12) = 2.0 * t * (r + s - 1.0);
		m(2, 13) = (-2.0) * r * t;
		m(2, 14) = (-2.0) * s * t;
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


