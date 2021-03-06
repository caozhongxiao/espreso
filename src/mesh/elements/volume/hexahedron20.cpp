
#include "mesh/elements/element.h"
#include "mesh/mesh.h"

#include "basis/containers/serializededata.h"
#include "basis/matrices/denseMatrix.h"

namespace espreso {

template<>
void Element::set<Element::CODE::HEXA20>()
{
	type = Element::TYPE::VOLUME;
	code = Element::CODE::HEXA20;
	nodes = 20;
	coarseNodes = 8;
	nCommonFace = 4;
	nCommonEdge = 3;

	std::vector<Element*> fpointers(6, &Mesh::edata[static_cast<int>(Element::CODE::SQUARE8)]);

	std::vector<int> data = {
		0, 1, 5, 4,  8, 17, 12, 16,
		3, 2, 1, 0, 10,  9,  8, 11,
		4, 5, 6, 7, 12, 13, 14, 15,
		7, 6, 2, 3, 14, 18, 10, 19,
		1, 2, 6, 5,  9, 18, 13, 17,
		3, 0, 4, 7, 11, 16, 15, 19
	};

	faces = new serializededata<int, int>(8, data);
	facepointers = new serializededata<int, Element*>(1, fpointers);

	size_t GPCount = 8, nodeCount = 20;

	std::vector<std::vector<double> > rst(3);

	switch (GPCount) {
	case 8: {
		double v = 0.577350269189625953;
		rst[0] = {  v,  v,  v,  v, -v, -v, -v, -v };
		rst[1] = { -v, -v,  v,  v, -v, -v,  v,  v };
		rst[2] = { -v,  v, -v,  v, -v,  v, -v,  v };
		break;
	}
	case 14: {
		double v1 = 0.758786910639329015;
		double v2 = 0.795822425754222018;
		double v3 = 0;
		rst[0] = { -v1,  v1,  v1, -v1, -v1,  v1,  v1, -v1,  v3,  v3,  v2, v3, -v2, v3 };
		rst[1] = { -v1, -v1,  v1,  v1, -v1, -v1,  v1,  v1,  v3, -v2,  v3, v2,  v3, v3 };
		rst[2] = { -v1, -v1, -v1, -v1,  v1,  v1,  v1,  v1, -v2,  v3,  v3, v3,  v3, v2 };
		break;
	}
	default:
		exit(1);
	}

	N = new std::vector<DenseMatrix>(GPCount, DenseMatrix(1, nodeCount));
	dN = new std::vector<DenseMatrix>(GPCount, DenseMatrix(3, nodeCount));
	weighFactor = new std::vector<double>();

	for (unsigned int i = 0; i < GPCount; i++) {
		DenseMatrix &m = (*N)[i];

		double r = rst[0][i];
		double s = rst[1][i];
		double t = rst[2][i];

		// basis function
		m(0, 0) = 0.125 * ((1.0 - r) * (1.0 - s) * (1.0 - t) * (-r - s - t - 2.0));
		m(0, 1) = 0.125 * ((1.0 + r) * (1.0 - s) * (1.0 - t) * ( r - s - t - 2.0));
		m(0, 2) = 0.125 * ((1.0 + r) * (1.0 + s) * (1.0 - t) * ( r + s - t - 2.0));
		m(0, 3) = 0.125 * ((1.0 - r) * (1.0 + s) * (1.0 - t) * (-r + s - t - 2.0));
		m(0, 4) = 0.125 * ((1.0 - r) * (1.0 - s) * (1.0 + t) * (-r - s + t - 2.0));
		m(0, 5) = 0.125 * ((1.0 + r) * (1.0 - s) * (1.0 + t) * ( r - s + t - 2.0));
		m(0, 6) = 0.125 * ((1.0 + r) * (1.0 + s) * (1.0 + t) * ( r + s + t - 2.0));
		m(0, 7) = 0.125 * ((1.0 - r) * (1.0 + s) * (1.0 + t) * (-r + s + t - 2.0));

		m(0, 8) =  0.25 * ((1.0 - pow(r, 2.0)) * (1.0 - s) * (1.0 - t));
		m(0, 9) =  0.25 * ((1.0 + r) * (1.0 - pow(s, 2.0)) * (1.0 - t));
		m(0, 10) = 0.25 * ((1.0 - pow(r, 2.0)) * (1.0 + s) * (1.0 - t));
		m(0, 11) = 0.25 * ((1.0 - r) * (1.0 - pow(s, 2.0)) * (1.0 - t));
		m(0, 12) = 0.25 * ((1.0 - pow(r, 2.0)) * (1.0 - s) * (1.0 + t));
		m(0, 13) = 0.25 * ((1.0 + r) * (1.0 - pow(s, 2.0)) * (1.0 + t));
		m(0, 14) = 0.25 * ((1.0 - pow(r, 2.0)) * (1.0 + s) * (1.0 + t));
		m(0, 15) = 0.25 * ((1.0 - r) * (1.0 - pow(s, 2.0)) * (1.0 + t));
		m(0, 16) = 0.25 * ((1.0 - r) * (1.0 - s) * (1.0 - pow(t, 2.0)));
		m(0, 17) = 0.25 * ((1.0 + r) * (1.0 - s) * (1.0 - pow(t, 2.0)));
		m(0, 18) = 0.25 * ((1.0 + r) * (1.0 + s) * (1.0 - pow(t, 2.0)));
		m(0, 19) = 0.25 * ((1.0 - r) * (1.0 + s) * (1.0 - pow(t, 2.0)));
	}

	for (unsigned int i = 0; i < GPCount; i++) {
		DenseMatrix &m = (*dN)[i];

		double r = rst[0][i];
		double s = rst[1][i];
		double t = rst[2][i];


		// dNr - derivation of basis function
		m(0, 0) =  ((s - 1.0) * (t - 1.0) * (r + s + t + 2.0)) / 8.0 + ((r - 1.0) * (s - 1.0) *         (t - 1.0)) / 8.0;
		m(0, 1) =  ((r + 1.0) * (s - 1.0) *         (t - 1.0)) / 8.0 - ((s - 1.0) * (t - 1.0) * (s - r + t + 2.0)) / 8.0;
		m(0, 2) = -((s + 1.0) * (t - 1.0) * (r + s - t - 2.0)) / 8.0 - ((r + 1.0) * (s + 1.0) *         (t - 1.0)) / 8.0;
		m(0, 3) = -((s + 1.0) * (t - 1.0) * (r - s + t + 2.0)) / 8.0 - ((r - 1.0) * (s + 1.0) *         (t - 1.0)) / 8.0;
		m(0, 4) = -((s - 1.0) * (t + 1.0) * (r + s - t + 2.0)) / 8.0 - ((r - 1.0) * (s - 1.0) *         (t + 1.0)) / 8.0;
		m(0, 5) = -((s - 1.0) * (t + 1.0) * (r - s + t - 2.0)) / 8.0 - ((r + 1.0) * (s - 1.0) *         (t + 1.0)) / 8.0;
		m(0, 6) =  ((s + 1.0) * (t + 1.0) * (r + s + t - 2.0)) / 8.0 + ((r + 1.0) * (s + 1.0) *         (t + 1.0)) / 8.0;
		m(0, 7) =  ((s + 1.0) * (t + 1.0) * (r - s - t + 2.0)) / 8.0 + ((r - 1.0) * (s + 1.0) *         (t + 1.0)) / 8.0;

		m(0, 8)  = -(r * (s - 1.0) * (t - 1.0)) / 2.0;
		m(0, 9)  =  ((pow(s, 2.0) - 1.0) * (t - 1.0)) / 4.0;
		m(0, 10) =  (r * (s + 1.0) * (t - 1.0)) / 2.0;
		m(0, 11) = -((pow(s, 2.0) - 1.0) * (t - 1.0)) / 4.0;
		m(0, 12) =  (r * (s - 1.0) * (t + 1.0)) / 2.0;
		m(0, 13) = -((pow(s, 2.0) - 1.0) * (t + 1.0)) / 4.0;
		m(0, 14) = -(r * (s + 1.0) * (t + 1.0)) / 2.0;
		m(0, 15) =  ((pow(s, 2.0) - 1.0) * (t + 1.0)) / 4.0;
		m(0, 16) = -((pow(t, 2.0) - 1.0) * (s - 1.0)) / 4.0;
		m(0, 17) =  ((pow(t, 2.0) - 1.0) * (s - 1.0)) / 4.0;
		m(0, 18) = -((pow(t, 2.0) - 1.0) * (s + 1.0)) / 4.0;
		m(0, 19) =  ((pow(t, 2.0) - 1.0) * (s + 1.0)) / 4.0;


		// dNs - derivation of basis function
		m(1, 0) =  ((r - 1.0) * (t - 1.0) * (r + s + t + 2.0)) / 8.0 + ((r - 1.0) * (s - 1.0) *         (t - 1.0)) / 8.0;
		m(1, 1) = -((r + 1.0) * (t - 1.0) * (s - r + t + 2.0)) / 8.0 - ((r + 1.0) * (s - 1.0) *         (t - 1.0)) / 8.0;
		m(1, 2) = -((r + 1.0) * (t - 1.0) * (r + s - t - 2.0)) / 8.0 - ((r + 1.0) * (s + 1.0) *         (t - 1.0)) / 8.0;
		m(1, 3) =  ((r - 1.0) * (s + 1.0) *         (t - 1.0)) / 8.0 - ((r - 1.0) * (t - 1.0) * (r - s + t + 2.0)) / 8.0;
		m(1, 4) = -((r - 1.0) * (t + 1.0) * (r + s - t + 2.0)) / 8.0 - ((r - 1.0) * (s - 1.0) *         (t + 1.0)) / 8.0;
		m(1, 5) =  ((r + 1.0) * (s - 1.0) *         (t + 1.0)) / 8.0 - ((r + 1.0) * (t + 1.0) * (r - s + t - 2.0)) / 8.0;
		m(1, 6) =  ((r + 1.0) * (t + 1.0) * (r + s + t - 2.0)) / 8.0 + ((r + 1.0) * (s + 1.0) *         (t + 1.0)) / 8.0;
		m(1, 7) =  ((r - 1.0) * (t + 1.0) * (r - s - t + 2.0)) / 8.0 - ((r - 1.0) * (s + 1.0) *         (t + 1.0)) / 8.0;

		m(1, 8)  = -((pow(r, 2.0) - 1.0) * (t - 1.0)) / 4.0;
		m(1, 9)  =  (s * (r + 1.0) * (t - 1.0)) / 2.0;
		m(1, 10) =  ((pow(r, 2.0) - 1.0) * (t - 1.0)) / 4.0;
		m(1, 11) = -(s * (r - 1.0) * (t - 1.0)) / 2.0;
		m(1, 12) =  ((pow(r, 2.0) - 1.0) * (t + 1.0)) / 4.0;
		m(1, 13) = -(s * (r + 1.0) * (t + 1.0)) / 2.0;
		m(1, 14) = -((pow(r, 2.0) - 1.0) * (t + 1.0)) / 4.0;
		m(1, 15) =  (s * (r - 1.0) * (t + 1.0)) / 2.0;
		m(1, 16) = -((pow(t, 2.0) - 1.0) * (r - 1.0)) / 4.0;
		m(1, 17) =  ((pow(t, 2.0) - 1.0) * (r + 1.0)) / 4.0;
		m(1, 18) = -((pow(t, 2.0) - 1.0) * (r + 1.0)) / 4.0;
		m(1, 19) =  ((pow(t, 2.0) - 1.0) * (r - 1.0)) / 4.0;

		// dNt - derivation of basis function
		m(2, 0) =  ((r - 1.0) * (s - 1.0) * (r + s + t + 2.0)) / 8.0 + ((r - 1.0) * (s - 1.0) *         (t - 1.0)) / 8.0;
		m(2, 1) = -((r + 1.0) * (s - 1.0) * (s - r + t + 2.0)) / 8.0 - ((r + 1.0) * (s - 1.0) *         (t - 1.0)) / 8.0;
		m(2, 2) =  ((r + 1.0) * (s + 1.0) *         (t - 1.0)) / 8.0 - ((r + 1.0) * (s + 1.0) * (r + s - t - 2.0)) / 8.0;
		m(2, 3) = -((r - 1.0) * (s + 1.0) * (r - s + t + 2.0)) / 8.0 - ((r - 1.0) * (s + 1.0) *         (t - 1.0)) / 8.0;
		m(2, 4) =  ((r - 1.0) * (s - 1.0) *         (t + 1.0)) / 8.0 - ((r - 1.0) * (s - 1.0) * (r + s - t + 2.0)) / 8.0;
		m(2, 5) = -((r + 1.0) * (s - 1.0) * (r - s + t - 2.0)) / 8.0 - ((r + 1.0) * (s - 1.0) *         (t + 1.0)) / 8.0;
		m(2, 6) =  ((r + 1.0) * (s + 1.0) * (r + s + t - 2.0)) / 8.0 + ((r + 1.0) * (s + 1.0) *         (t + 1.0)) / 8.0;
		m(2, 7) =  ((r - 1.0) * (s + 1.0) * (r - s - t + 2.0)) / 8.0 - ((r - 1.0) * (s + 1.0) *         (t + 1.0)) / 8.0;

		m(2, 8)  = -((pow(r, 2.0) - 1.0) * (s - 1.0)) / 4.0;
		m(2, 9)  =  ((pow(s, 2.0) - 1.0) * (r + 1.0)) / 4.0;
		m(2, 10) =  ((pow(r, 2.0) - 1.0) * (s + 1.0)) / 4.0;
		m(2, 11) = -((pow(s, 2.0) - 1.0) * (r - 1.0)) / 4.0;
		m(2, 12) =  ((pow(r, 2.0) - 1.0) * (s - 1.0)) / 4.0;
		m(2, 13) = -((pow(s, 2.0) - 1.0) * (r + 1.0)) / 4.0;
		m(2, 14) = -((pow(r, 2.0) - 1.0) * (s + 1.0)) / 4.0;
		m(2, 15) =  ((pow(s, 2.0) - 1.0) * (r - 1.0)) / 4.0;
		m(2, 16) = -(t * (r - 1.0) * (s - 1.0)) / 2.0;
		m(2, 17) =  (t * (r + 1.0) * (s - 1.0)) / 2.0;
		m(2, 18) = -(t * (r + 1.0) * (s + 1.0)) / 2.0;
		m(2, 19) =  (t * (r - 1.0) * (s + 1.0)) / 2.0;
	}

	switch (GPCount) {
	case 8: {
		weighFactor->resize(8, 1.0);
		break;
	}
	case 14: {
		weighFactor->resize(8, 0.335180055401662);
		weighFactor->resize(14, 0.886426592797784);
		break;
	}
	default:
		exit(1);
	}
}
}



