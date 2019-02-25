
#include "mesh/elements/element.h"
#include "mesh/mesh.h"

#include "basis/containers/serializededata.h"
#include "basis/matrices/denseMatrix.h"

namespace espreso {

template<>
void Element::set<Element::CODE::PYRAMID5>()
{
	type = Element::TYPE::VOLUME;
	code = Element::CODE::PYRAMID5;
	nodes = 5;
	coarseNodes = 5;
	nCommonFace = 3;
	nCommonEdge = 2;

	std::vector<Element*> fpointers;
	fpointers.resize(1, &Mesh::edata[static_cast<int>(Element::CODE::SQUARE4)]);
	fpointers.resize(5, &Mesh::edata[static_cast<int>(Element::CODE::TRIANGLE3)]);

	std::vector<int> dist = { 0, 4, 7, 10, 13, 16 };
	std::vector<int> data = {
		3, 2, 1, 0,
		0, 1, 4,
		1, 2, 4,
		2, 3, 4,
		3, 0, 4
	};

	faces = new serializededata<int, int>(dist, data);
	facepointers = new serializededata<int, Element*>(1, fpointers);

	size_t GPCount = 8, nodeCount = 5;

	N = new std::vector<DenseMatrix>(GPCount, DenseMatrix(1, nodeCount));
	dN = new std::vector<DenseMatrix>(GPCount, DenseMatrix(3, nodeCount));
	weighFactor = new std::vector<double>();

	std::vector< std::vector<double> > rst(3, std::vector<double>(GPCount));

	switch (GPCount) {
		case 8: {
			double v = 0.577350269189625953;
			rst[0] = {  v,  v,  v,  v, -v, -v, -v, -v };
			rst[1] = { -v, -v,  v,  v, -v, -v,  v,  v };
			rst[2] = { -v,  v, -v,  v, -v,  v, -v,  v };
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
		m(0, 0) = 0.125 * ((1 - r) * (1 - s) * (1 - t));
		m(0, 1) = 0.125 * ((1 + r) * (1 - s) * (1 - t));
		m(0, 2) = 0.125 * ((1 + r) * (1 + s) * (1 - t));
		m(0, 3) = 0.125 * ((1 - r) * (1 + s) * (1 - t));
		m(0, 4) = 0.125 * ( 4 * (1 + t));
	}

	for (unsigned int i = 0; i < GPCount; i++) {
		///dN contains [dNr, dNs, dNt]
		DenseMatrix &m = (*dN)[i];

		double r = rst[0][i];
		double s = rst[1][i];
		double t = rst[2][i];

		// dNr - derivation of basis function
		m(0, 0) = 0.125 * (-(1. - s) * (1. - t));
		m(0, 1) = 0.125 * ( (1. - s) * (1. - t));
		m(0, 2) = 0.125 * ( (1. + s) * (1. - t));
		m(0, 3) = 0.125 * (-(1. + s) * (1. - t));
		m(0, 4) = 0;

		// dNs - derivation of basis function
		m(1, 0) = 0.125 * (-(1. - r) * (1. - t));
		m(1, 1) = 0.125 * (-(1. + r) * (1. - t));
		m(1, 2) = 0.125 * ( (1. + r) * (1. - t));
		m(1, 3) = 0.125 * ( (1. - r) * (1. - t));
		m(1, 4) = 0;

		// dNt - derivation of basis function
		m(2, 0) = 0.125 * (-(1. - r) * (1. - s));
		m(2, 1) = 0.125 * (-(1. + r) * (1. - s));
		m(2, 2) = 0.125 * (-(1. + r) * (1. + s));
		m(2, 3) = 0.125 * (-(1. - r) * (1. + s));
		m(2, 4) = 0.125 * (4.0);
	}

	switch (GPCount) {
	case 8: {
		(*weighFactor) = std::vector<double> (8, 1.0);
		break;
	}
	default:
		exit(1);
	}
}
}


