
#include "tetrahedron10.h"

#include "basis/containers/serializededata.h"
#include "basis/matrices/denseMatrix.h"

using namespace espreso;

Element Tetrahedron10::fill(Element e, Element* begin)
{
	std::vector<Element*> facepointers(4, begin + static_cast<int>(Element::CODE::TRIANGLE6));

	std::vector<int> data = {
		0, 1, 3, 4, 8, 7,
		1, 2, 3, 5, 9, 8,
		2, 0, 3, 6, 7, 9,
		2, 1, 0, 5, 4, 6
	};

	e.faces = new serializededata<int, int>(6, data);
	e.facepointers = new serializededata<int, Element*>(1, facepointers);

	size_t GPCount = 15, nodeCount = 10;

	std::vector<std::vector<double> > rst(3);

	switch (GPCount) {
	case 15: {
		rst[0] = {
				0.2500000000000000, 0.0000000000000000, 0.3333333333333333, 0.3333333333333333,
				0.3333333333333333, 0.7272727272727273, 0.0909090909090909, 0.0909090909090909,
				0.0909090909090909, 0.4334498464263357, 0.0665501535736643, 0.0665501535736643,
				0.0665501535736643, 0.4334498464263357, 0.4334498464263357};
		rst[1] = {
				0.2500000000000000, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
				0.0000000000000000, 0.0909090909090909, 0.0909090909090909, 0.0909090909090909,
				0.7272727272727273, 0.0665501535736643, 0.4334498464263357, 0.0665501535736643,
				0.4334498464263357, 0.0665501535736643, 0.4334498464263357};
		rst[2] = {
				0.2500000000000000, 0.3333333333333333, 0.3333333333333333, 0.0000000000000000,
				0.3333333333333333, 0.0909090909090909, 0.0909090909090909, 0.7272727272727273,
				0.0909090909090909, 0.0665501535736643, 0.0665501535736643, 0.4334498464263357,
				0.4334498464263357, 0.4334498464263357, 0.0665501535736643};
		break;
	}
	default:
		exit(1);
	}

	e.N = new std::vector<DenseMatrix>(GPCount, DenseMatrix(1, nodeCount));
	e.dN = new std::vector<DenseMatrix>(GPCount, DenseMatrix(3, nodeCount));
	e.weighFactor = new std::vector<double>();

	for (unsigned int i = 0; i < GPCount; i++) {
		double r = rst[0][i];
		double s = rst[1][i];
		double t = rst[2][i];

		DenseMatrix &m = (*e.N)[i];

		m(0, 0) = r * (2.0 * r - 1.0);
		m(0, 1) = s * (2.0 * s - 1.0);
		m(0, 2) = t * (2.0 * t - 1.0);
		m(0, 3) = 2.0 * r * r + 4.0 * r * s + 4.0 * r * t - 3.0 * r + 2.0* s * s + 4.0 * s * t - 3.0 * s + 2.0 * t * t - 3.0 * t + 1.0;
		m(0, 4) = 4.0 * r * s;
		m(0, 5) = 4.0 * s * t;
		m(0, 6) = 4.0 * r * t;
		m(0, 7) = r * (-4.0 * r - 4.0 * s - 4.0 * t + 4.0);
		m(0, 8) = s * (-4.0 * r - 4.0 * s - 4.0 * t + 4.0);
		m(0, 9) = t * (-4.0 * r - 4.0 * s - 4.0 * t + 4.0);
	}

	for (unsigned int i = 0; i < GPCount; i++) {
		double r = rst[0][i];
		double s = rst[1][i];
		double t = rst[2][i];

		DenseMatrix &m = (*e.dN)[i];

		m(0, 0) = 4.0 * r - 1.0;
		m(0, 1) = 0;
		m(0, 2) = 0;
		m(0, 3) = 4.0 * r + 4.0 * s + 4.0 * t - 3.0;
		m(0, 4) = 4.0 * s;
		m(0, 5) = 0;
		m(0, 6) = 4.0 * t;
		m(0, 7) = -8.0 * r - 4.0 * s - 4.0 * t + 4.0;
		m(0, 8) = -4.0 * s;
		m(0, 9) = -4.0 * t;

		m(2, 0) = 0;
		m(2, 1) = 4.0 * s - 1.0;
		m(2, 2) = 0 ;
		m(2, 3) = 4.0 * r + 4.0 * s + 4.0 * t - 3.0;
		m(2, 4) = 4.0 * r;
		m(2, 5) = 4.0 * t;
		m(2, 6) = 0;
		m(2, 7) = -4.0 * r;
		m(2, 8) = -4.0 * r - 8.0 * s - 4.0 * t + 4.0;
		m(2, 9) = -4.0 * t;

		m(1, 0) = 0;
		m(1, 1) = 0;
		m(1, 2) = 4.0 * t - 1.0;
		m(1, 3) = 4.0 * r + 4.0 * s + 4.0* t  - 3.0;
		m(1, 4) = 0;
		m(1, 5) = 4.0 * s;
		m(1, 6) = 4.0 * r;
		m(1, 7) = -4.0 * r;
		m(1, 8) = -4.0 * s;
		m(1, 9) = -4.0 * r - 4.0 * s - 8.0 * t + 4.0;
	}

	switch (GPCount) {
	case 15: {
		e.weighFactor->resize( 1, 0.030283678097089);
		e.weighFactor->resize( 5, 0.006026785714286);
		e.weighFactor->resize( 9, 0.011645249086029);
		e.weighFactor->resize(15, 0.010949141561386);
		break;
	}
	default:
		exit(1);
	}

	return e;
}







