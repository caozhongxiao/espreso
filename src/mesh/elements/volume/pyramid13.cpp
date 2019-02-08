
#include "mesh/elements/element.h"
#include "mesh/mesh.h"

#include "basis/containers/serializededata.h"
#include "basis/matrices/denseMatrix.h"

using namespace espreso;

template<>
void Element::set<Element::CODE::PYRAMID13>()
{
	type = Element::TYPE::VOLUME;
	code = Element::CODE::PYRAMID13;
	nodes = 13;
	coarseNodes = 5;
	nCommonFace = 4;
	nCommonEdge = 3;

	std::vector<Element*> fpointers;
	fpointers.resize(1, &Mesh::edata[static_cast<int>(Element::CODE::SQUARE8)]);
	fpointers.resize(5, &Mesh::edata[static_cast<int>(Element::CODE::TRIANGLE6)]);

	std::vector<int> dist = { 0, 8, 14, 20, 26, 32 };
	std::vector<int> data = {
		3, 2, 1, 0,  7,  6, 5, 8,
		0, 1, 4, 5, 10,  9,
		1, 2, 4, 6, 11, 10,
		2, 3, 4, 7, 12, 11,
		3, 0, 4, 8,  9, 12
	};

	faces = new serializededata<int, int>(dist, data);
	facepointers = new serializededata<int, Element*>(1, fpointers);

	size_t GPCount = 14, nodeCount = 13;

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

	for (unsigned int i = 0; i < GPCount; i++) {
		DenseMatrix &m = (*N)[i];

		double r = rst[0][i];
		double s = rst[1][i];
		double t = rst[2][i];

		// basis function
		m(0, 0)  = ((0.5 * (1.0 - t)) / 4.0) * ((1.0 - r) * (1.0 - s) * (-1.0 - (0.5 * (1.0 - t)) * r - (0.5 * (1.0 - t)) * s));
		m(0, 1)  = ((0.5 * (1.0 - t)) / 4.0) * ((1.0 + r) * (1.0 - s) * (-1.0 + (0.5 * (1 - t))   * r - (0.5 * (1.0 - t)) * s));
		m(0, 2)  = ((0.5 * (1.0 - t)) / 4.0) * ((1.0 + r) * (1.0 + s) * (-1.0 + (0.5 * (1.0 - t)) * r + (0.5 * (1.0 - t)) * s));
		m(0, 3)  = ((0.5 * (1.0 - t)) / 4.0) * ((1.0 - r) * (1.0 + s) * (-1.0 - (0.5 * (1.0 - t)) * r + (0.5 * (1.0 - t)) * s));
		m(0, 4)  = (1.0 - (0.5 * (1.0 - t))) * (1.0 - 2.0 * (0.5 * (1.0 - t)));
		m(0, 5)  = (((0.5 * (1.0 - t)) * (0.5 * (1.0 - t))) / 2.0) * (1.0 - s) * (1.0 - pow(r, 2.0));
		m(0, 6)  = (((0.5 * (1.0 - t)) * (0.5 * (1.0 - t))) / 2.0) * (1.0 + r) * (1.0 - pow(s, 2.0));
		m(0, 7)  = (((0.5 * (1.0 - t)) * (0.5 * (1.0 - t))) / 2.0) * (1.0 + s) * (1.0 - pow(r, 2.0));
		m(0, 8)  = (((0.5 * (1.0 - t)) * (0.5 * (1.0 - t))) / 2.0) * (1.0 - r) * (1.0 - pow(s, 2.0));
		m(0, 9)  = (0.5 * (1.0 - t)) * (1.0 - (0.5 * (1.0 - t))) * (1.0 - r - s + r * s);
		m(0, 10) = (0.5 * (1.0 - t)) * (1.0 - (0.5 * (1.0 - t))) * (1.0 + r - s - r * s);
		m(0, 11) = (0.5 * (1.0 - t)) * (1.0 - (0.5 * (1.0 - t))) * (1.0 + r + s + r * s);
		m(0, 12) = (0.5 * (1.0 - t)) * (1.0 - (0.5 * (1.0 - t))) * (1.0 - r + s - r * s);
	}

	for (unsigned int i = 0; i < GPCount; i++) {
		///dN contains [dNr, dNs, dNt]
		DenseMatrix &m = (*dN)[i];

		double r = rst[0][i];
		double s = rst[1][i];
		double t = rst[2][i];

		// dNr - derivation of basis function
		m(0, 0)  = -(t / 8.0 - 1.0 / 8.0) * (s - 1.0) * (r * (t / 2.0 - 0.5) + s * (t / 2.0 - 0.5) - 1.0) - (t / 2.0 - 0.5) * (t / 8.0 - 1.0 / 8.0) * (r - 1.0) * (s - 1.0);
		m(0, 1)  = -(t / 8.0 - 1.0 / 8.0) * (s - 1.0) * (r * (t / 2.0 - 0.5) - s * (t / 2.0 - 0.5) + 1.0) - (t / 2.0 - 0.5) * (t / 8.0 - 1.0 / 8.0) * (r + 1.0) * (s - 1.0);
		m(0, 2)  =  (t / 8.0 - 1.0 / 8.0) * (s + 1.0) * (r * (t / 2.0 - 0.5) + s * (t / 2.0 - 0.5) + 1.0) + (t / 2.0 - 0.5) * (t / 8.0 - 1.0 / 8.0) * (r + 1.0) * (s + 1.0);
		m(0, 3)  =  (t / 2.0 - 0.5) * (s + 1.0) * (r - 1.0) * (t / 8.0 - 1.0 / 8.0) - (t / 8.0 - 1.0 / 8.0) * (s + 1.0) * (s * (t / 2.0 - 0.5) - r * (t / 2.0 - 0.5) + 1.0);
		m(0, 4)  =  0.0;
		m(0, 5)  =  r * ((t / 2.0 - 0.5) * (t / 2.0 - 0.5)) * (s - 1.0);
		m(0, 6)  = -((pow(s, 2.0) - 1.0) * (t / 2.0 - 0.5) * (t / 2.0 - 0.5)) / 2.0;
		m(0, 7)  = -r * ((t / 2.0 - 0.5) * (t / 2.0 - 0.5)) * (s + 1.0);
		m(0, 8)  =  ((pow(s, 2.0) - 1) * (t / 2.0 - 0.5) * (t / 2.0 - 0.5)) / 2.0;
		m(0, 9)  = -(t / 2.0 - 0.5) * (t / 2.0 + 0.5) * (s - 1.0);
		m(0, 10) =  (t / 2.0 - 0.5) * (t / 2.0 + 0.5) * (s - 1.0);
		m(0, 11) = -(t / 2.0 - 0.5) * (t / 2.0 + 0.5) * (s + 1.0);
		m(0, 12) =  (t / 2.0 - 0.5) * (t / 2.0 + 0.5) * (s + 1.0);
		//  dNs - derivation of basis function
		m(1, 0)  = -(t / 8.0 - 1.0 / 8.0) * (r - 1.0) * (r * (t / 2.0 - 1.0 / 2.0) + s * (t / 2.0 - 1.0 / 2.0) - 1.0) - (t / 2.0 - 1.0 / 2.0) * (t / 8.0 - 1.0 / 8.0) * (r - 1.0) * (s - 1.0);
		m(1, 1)  =  (t / 2.0 - 1.0 / 2.0) * (t / 8.0 - 1.0 / 8.0) * (r + 1.0) * (s - 1.0) - (t / 8.0 - 1.0 / 8.0) * (r + 1.0) * (r * (t / 2.0 - 1.0 / 2.0) - s * (t / 2.0 - 1.0 / 2.0) + 1.0);
		m(1, 2)  =  (t / 8.0 - 1.0 / 8.0) * (r + 1.0) * (r * (t / 2.0 - 1.0 / 2.0) + s * (t / 2.0 - 1.0 / 2.0) + 1.0) + (t / 2.0 - 1.0 / 2.0) * (t / 8.0 - 1.0 / 8.0) * (r + 1.0) * (s + 1.0);
		m(1, 3)  = -(t / 8.0 - 1.0 / 8.0) * (r - 1.0) * (s * (t / 2.0 - 1.0 / 2.0) - r * (t / 2.0 - 1.0 / 2.0) + 1.0) - (t / 2.0 - 1.0 / 2.0) * (t / 8.0 - 1.0 / 8.0) * (r - 1.0) * (s + 1.0);
		m(1, 4)  =  0.0;
		m(1, 5)  =  ((pow(r, 2.0) - 1.0) * (t / 2.0 - 0.5) * (t / 2.0 - 1.0 / 2.0)) / 2.0;
		m(1, 6)  = -s * (t / 2.0 - 1.0 / 2.0) * (t / 2.0 - 0.5) * (r + 1.0);
		m(1, 7)  = -((pow(r, 2.0) - 1.0) * (t / 2.0 - 0.5) * (t / 2.0 - 1.0 / 2.0)) / 2.0;
		m(1, 8)  =  s * (t / 2.0 - 0.5) * (t / 2.0 - 0.5) * (r - 1.0);
		m(1, 9)  = -(t / 2.0 - 0.5) * (t / 2.0 + 0.5) * (r - 1.0);
		m(1, 10) =  (t / 2.0 - 0.5) * (t / 2.0 + 0.5) * (r + 1.0);
		m(1, 11) = -(t / 2.0 - 0.5) * (t / 2.0 + 0.5) * (r + 1.0);
		m(1, 12) =  (t / 2.0 - 0.5) * (t / 2.0 + 0.5) * (r - 1.0);
		//  dNt - derivation of basis function
		m(2, 0)  = -((r - 1.0) * (s - 1.0) * (r * (t / 2.0 - 0.5) + s * (t / 2.0 - 0.5) - 1.0)) / 8.0 - (r / 2.0 + s / 2.0) * (t / 8.0 - 1.0 / 8.0) * (r - 1.0) * (s - 1.0);
		m(2, 1)  = -((r + 1.0) * (s - 1.0) * (r * (t / 2.0 - 0.5) - s * (t / 2.0 - 0.5) + 1.0)) / 8.0 - (r / 2.0 - s / 2.0) * (t / 8.0 - 1.0 / 8.0) * (r + 1.0) * (s - 1.0);
		m(2, 2)  =  ((r + 1.0) * (s + 1.0) * (r * (t / 2.0 - 0.5) + s * (t / 2.0 - 0.5) + 1.0)) / 8.0 + (r / 2.0 + s / 2.0) * (t / 8.0 - 1.0 / 8.0) * (r + 1.0) * (s + 1.0);
		m(2, 3)  =  (r / 2.0 - s / 2.0) * (t / 8.0 - 1.0 / 8.0) * (r - 1.0) * (s + 1.0) - ((r - 1.0) * (s + 1.0) * (s * (t / 2.0 - 0.5) - r * (t / 2.0 - 0.5) + 1.0)) / 8.0;
		m(2, 4)  =  t + 0.5;
		m(2, 5)  =  ((pow(r, 2.0) - 1.0) * (t / 2.0 - 0.5) * (s - 1.0)) / 2.0;
		m(2, 6)  = -((pow(s, 2.0) - 1.0) * (t / 2.0 - 0.5) * (r + 1.0)) / 2.0;
		m(2, 7)  = -((pow(r, 2.0) - 1.0) * (t / 2.0 - 0.5) * (s + 1.0)) / 2.0;
		m(2, 8)  =  ((pow(s, 2.0) - 1.0) * (t / 2.0 - 0.5) * (r - 1.0)) / 2.0;
		m(2, 9)  =  ((t / 2.0 - 0.5) * (r + s - r * s - 1.0)) / 2.0 + ((t / 2.0 + 0.5) * (r + s - r * s - 1.0)) / 2.0;
		m(2, 10) = -((t / 2.0 - 0.5) * (r - s - r * s + 1.0)) / 2.0 - ((t / 2.0 + 0.5) * (r - s - r * s + 1.0)) / 2.0;
		m(2, 11) = -((t / 2.0 - 0.5) * (r + s + r * s + 1.0)) / 2.0 - ((t / 2.0 + 0.5) * (r + s + r * s + 1.0)) / 2.0;
		m(2, 12) =  ((t / 2.0 - 0.5) * (r - s + r * s - 1.0)) / 2.0 + ((t / 2.0 + 0.5) * (r - s + r * s - 1.0)) / 2.0;
	}

	switch (GPCount) {
	case 8: {
		(*weighFactor) = std::vector<double> (8, 1.0);
		break;
	}
	case 14: {
		(*weighFactor).resize(8, 0.335180055401662);
		(*weighFactor).resize(14, 0.886426592797784);
		break;
	}
	default:
		exit(1);
	}
}







