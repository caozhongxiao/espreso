
#include "mesh/elements/element.h"
#include "mesh/mesh.h"

#include "basis/containers/serializededata.h"
#include "basis/matrices/denseMatrix.h"

namespace espreso {

template<>
void Element::set<Element::CODE::TETRA4>()
{
	type = Element::TYPE::VOLUME;
	code = Element::CODE::TETRA4;
	nodes = 4;
	coarseNodes = 4;
	nCommonFace = 3;
	nCommonEdge = 2;

	std::vector<Element*> fpointers(4, &Mesh::edata[static_cast<int>(Element::CODE::TRIANGLE3)]);

	std::vector<int> data = {
		0, 1, 3,
		1, 2, 3,
		2, 0, 3,
		2, 1, 0
	};

	faces = new serializededata<int, int>(3, data);
	facepointers = new serializededata<int, Element*>(1, fpointers);

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

	N = new std::vector<DenseMatrix>(GPCount, DenseMatrix(1, nodeCount));
	dN = new std::vector<DenseMatrix>(GPCount, DenseMatrix(3, nodeCount));
	weighFactor = new std::vector<double>();

	for (unsigned int i = 0; i < GPCount; i++) {
		double r = rst[0][i];
		double s = rst[1][i];
		double t = rst[2][i];

		DenseMatrix &m = (*N)[i];

		m(0, 0) = r;
		m(0, 1) = s;
		m(0, 2) = t;
		m(0, 3) = 1.0 - r - s - t;
	}

	for (unsigned int i = 0; i < GPCount; i++) {
		//  N = [ r, s, t,  1 - r - s - t ];
		DenseMatrix &m = (*dN)[i];

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
		weighFactor->resize(4, 1.0 / 24.0);
		break;
	}
	default:
		exit(1);
	}
}
}






