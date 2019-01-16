
#include "tetrahedron4.h"

#include "basis/containers/serializededata.h"
#include "basis/matrices/denseMatrix.h"

using namespace espreso;

Element Tetrahedron4::fill(Element e, Element* begin)
{
	std::vector<Element*> facepointers(4, begin + static_cast<int>(Element::CODE::TRIANGLE3));

	std::vector<int> data = {
		0, 1, 3,
		1, 2, 3,
		2, 0, 3,
		2, 1, 0
	};

	e.faces = new serializededata<int, int>(3, data);
	e.facepointers = new serializededata<int, Element*>(1, facepointers);

	size_t GPCount = 4, nodeCount = 4;

	std::vector<std::vector<double> > rst(3);

	switch (GPCount) {
	case 4: {
		rst[0] = {0.5854101966249685, 0.1381966011250105, 0.1381966011250105, 0.1381966011250105};
		rst[1] = {0.1381966011250105, 0.1381966011250105, 0.1381966011250105, 0.5854101966249685};
		rst[2] = {0.1381966011250105, 0.1381966011250105, 0.5854101966249685, 0.1381966011250105};
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

		m(0, 0) = r;
		m(0, 1) = s;
		m(0, 2) = t;
		m(0, 3) = 1.0 - r - s - t;
	}

	for (unsigned int i = 0; i < GPCount; i++) {
		//  N = [ r, s, t,  1 - r - s - t ];
		DenseMatrix &m = (*e.dN)[i];

		// dNr = [ 1, 0, 0, -1 ];
		m(0, 0) =  1.0;
		m(0, 1) =  0.0;
		m(0, 2) =  0.0;
		m(0, 3) = -1.0;

		// dNs = [ 0, 1, 0, -1 ];
		m(2, 0) =  0.0;
		m(2, 1) =  1.0;
		m(2, 2) =  0.0;
		m(2, 3) = -1.0;

		// dNs = [ 0, 0, 1, -1 ];
		m(1, 0) =  0.0;
		m(1, 1) =  0.0;
		m(1, 2) =  1.0;
		m(1, 3) = -1.0;
	}

	switch (GPCount) {
	case 4: {
		e.weighFactor->resize(4, 1.0 / 24.0);
		break;
	}
	default:
		exit(1);
	}

	return e;
}







