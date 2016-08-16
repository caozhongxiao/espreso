
#include "tetrahedron4.h"
#include "../line/line2.h"
#include "../plane/triangle3.h"

using namespace espreso;

size_t Tetrahedron4::_counter = 0;

std::vector<Property> Tetrahedron4::_DOFElement;
std::vector<Property> Tetrahedron4::_DOFFace;
std::vector<Property> Tetrahedron4::_DOFEdge;
std::vector<Property> Tetrahedron4::_DOFPoint;
std::vector<Property> Tetrahedron4::_DOFMidPoint;

static std::vector<DenseMatrix> Tetra4_dN()
{
	// dN contains [dNr, dNs, dNt]
	std::vector<DenseMatrix> dN(
		Tetrahedron4GPCount,
		DenseMatrix(Point::size(), Tetrahedron4NodesCount)
	);

	for (unsigned int i = 0; i < Tetrahedron4GPCount; i++) {
		//  N = [ r, s, t,  1 - r - s - t ];
		DenseMatrix &m = dN[i];

		// dNr = [ 1, 0, 0, -1 ];
		m(0, 0) =  1.0;
		m(0, 1) =  0.0;
		m(0, 2) =  0.0;
		m(0, 3) = -1.0;

		// dNs = [ 0, 1, 0, -1 ];
		m(1, 0) =  0.0;
		m(1, 1) =  1.0;
		m(1, 2) =  0.0;
		m(1, 3) = -1.0;

		// dNs = [ 0, 0, 1, -1 ];
		m(2, 0) =  0.0;
		m(2, 1) =  0.0;
		m(2, 2) =  1.0;
		m(2, 3) = -1.0;
	}

	return dN;
}

static std::vector<DenseMatrix> Tetra4_N()
{
	std::vector<DenseMatrix> N(
			Tetrahedron4GPCount,
			DenseMatrix(1, Tetrahedron4NodesCount));

	std::vector<double> rv;
	std::vector<double> sv;
	std::vector<double> tv;

	switch (Tetrahedron4GPCount) {
	case 4:
		rv = {0.5854101966249685, 0.1381966011250105, 0.1381966011250105, 0.1381966011250105};
		sv = {0.1381966011250105, 0.1381966011250105, 0.1381966011250105, 0.5854101966249685};
		tv = {0.1381966011250105, 0.1381966011250105, 0.5854101966249685, 0.1381966011250105};
		break;
	case 5:
		rv = {0.2500000000000000, 0.5000000000000000, 0.1666666666666667, 0.1666666666666667, 0.1666666666666667};
		sv = {0.2500000000000000, 0.1666666666666667, 0.1666666666666667, 0.1666666666666667, 0.5000000000000000};
		tv = {0.2500000000000000, 0.1666666666666667, 0.1666666666666667, 0.5000000000000000, 0.1666666666666667};
		break;
	case 11:
		rv = {
				0.2500000000000000, 0.7857142857142857, 0.0714285714285714, 0.0714285714285714,
				0.0714285714285714, 0.1005964238332008, 0.3994035761667992, 0.3994035761667992,
				0.3994035761667992, 0.1005964238332008, 0.1005964238332008};
		sv = {
				0.2500000000000000, 0.0714285714285714, 0.0714285714285714, 0.0714285714285714,
				0.7857142857142857, 0.3994035761667992, 0.1005964238332008, 0.3994035761667992,
				0.1005964238332008, 0.3994035761667992, 0.1005964238332008};
		tv = {
				0.2500000000000000, 0.0714285714285714, 0.0714285714285714, 0.7857142857142857,
				0.0714285714285714, 0.3994035761667992, 0.3994035761667992, 0.1005964238332008,
				0.1005964238332008, 0.1005964238332008, 0.3994035761667992};
		break;
	case 15:
		rv = {
				0.2500000000000000, 0.0000000000000000, 0.3333333333333333, 0.3333333333333333,
				0.3333333333333333, 0.7272727272727273, 0.0909090909090909, 0.0909090909090909,
				0.0909090909090909, 0.4334498464263357, 0.0665501535736643, 0.0665501535736643,
				0.0665501535736643, 0.4334498464263357, 0.4334498464263357};
		sv = {
				0.2500000000000000, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
				0.0000000000000000, 0.0909090909090909, 0.0909090909090909, 0.0909090909090909,
				0.7272727272727273, 0.0665501535736643, 0.4334498464263357, 0.0665501535736643,
				0.4334498464263357, 0.0665501535736643, 0.4334498464263357};
		tv = {
				0.2500000000000000, 0.3333333333333333, 0.3333333333333333, 0.0000000000000000,
				0.3333333333333333, 0.0909090909090909, 0.0909090909090909, 0.7272727272727273,
				0.0909090909090909, 0.0665501535736643, 0.0665501535736643, 0.4334498464263357,
				0.4334498464263357, 0.4334498464263357, 0.0665501535736643};
		break;
	default:
		ESINFO(ERROR) << "Unknown number of Tatrahedron4 GP count.";
	}


	// N = [r s t  1-r-s-t];
	for (unsigned int i = 0; i < Tetrahedron4GPCount; i++) {
		double r = rv[i];
		double s = sv[i];
		double t = tv[i];

		N[i](0, 0) = r;
		N[i](0, 1) = s;
		N[i](0, 2) = t;
		N[i](0, 3) = 1.0 - r - s - t;
	}

	return N;
}

static std::vector<double> Tetra4_Weight()
{
	switch (Tetrahedron4GPCount) {
	case 4: {
		return std::vector<double> (4, 1.0 / 24.0);
	}
	case 5: {
		std::vector<double> w(5, 3.0 / 40.0);
		w[0] = - 2.0 / 15.0;
		return w;
	}
	case 11: {
		std::vector<double> w(11);
		w[0] = -0.013155555555556;
		w[1] = w[2] = w[3] = w[4] = 0.007622222222222;
		w[5] = w[6] = w[7] = w[8] = w[9] = w[10] = 0.024888888888889;
		return w;
	}
	case 15: {
		std::vector<double> w(15);
		w[0] = 0.030283678097089;
		w[1] = w[2] = w[3] = w[4] = 0.006026785714286;
		w[5] = w[6] = w[7] = w[8] = 0.011645249086029;
		w[9] = w[10] = w[11] = w[12] = w[13] = w[14] = 0.010949141561386;
		return w;
	}
	default:
		ESINFO(ERROR) << "Unknown number of Tatrahedron4 GP count.";
		exit(EXIT_FAILURE);
	}
}

std::vector<DenseMatrix> Tetrahedron4::_dN = Tetra4_dN();
std::vector<DenseMatrix> Tetrahedron4::_N = Tetra4_N();
std::vector<double> Tetrahedron4::_weighFactor = Tetra4_Weight();

bool Tetrahedron4::match(const eslocal *indices, eslocal n) {

#if ESPRESO_POINT_DIMENSION == 2
	// Tetrahedron4 is 3D element
	return false;
#endif

	switch (n) {
	case 4:
		for (eslocal i = 0; i < 3; i++) {
			for (eslocal j = i + 1; j < 4; j++) {
				if (Element::match(indices, i, j)) {
					return false;
				}
			}
		}
		return true;
	case 8:
		if (!Element::match(indices, 2, 3)) {
			return false;
		}
		if (!Element::match(indices, 4, 5)) {
			return false;
		}
		if (!Element::match(indices, 5, 6)) {
			return false;
		}
		if (!Element::match(indices, 6, 7)) {
			return false;
		}

		eslocal various[4] = { 0, 1, 2, 4 };
		for (eslocal i = 0; i < 3; i++) {
			for (eslocal j = i + 1; j < 4; j++) {
				if (Element::match(indices, various[i], various[j])) {
					return false;
				}
			}
		}
		return true;
	default:
		return false;
	}
}

std::vector<eslocal> Tetrahedron4::getNeighbours(size_t nodeIndex) const
{
	std::vector<eslocal> result;
	result.reserve(3);
	for (size_t i = 0; i < Tetrahedron4NodesCount; i++) {
		if (i != nodeIndex) {
			result.push_back(_indices[i]);
		}
	}
	return result;
}

void Tetrahedron4::fillEdges()
{
	eslocal line[Line2NodesCount];
	_edges.reserve(Tetrahedron4EdgeCount);

	for (size_t edge = 0; edge < 3; edge++) {
		line[0] = _indices[ edge         ];
		line[1] = _indices[(edge + 1) % 3];
		_edges.push_back(new Line2(line));
		_edges.back()->elements().push_back(this);

		line[0] = _indices[edge];
		line[1] = _indices[   3];
		_edges.push_back(new Line2(line));
		_edges.back()->elements().push_back(this);
	}
}

void Tetrahedron4::fillFaces()
{
	eslocal triangle[Triangle3NodesCount];
	_faces.reserve(Tetrahedron4FacesCount);

	triangle[0] = _indices[1];
	triangle[1] = _indices[0];
	triangle[2] = _indices[2];
	_faces.push_back(new Triangle3(triangle));
	_faces.back()->elements().push_back(this);

	triangle[0] = _indices[0];
	triangle[1] = _indices[1];
	triangle[2] = _indices[3];
	_faces.push_back(new Triangle3(triangle));
	_faces.back()->elements().push_back(this);

	triangle[0] = _indices[1];
	triangle[1] = _indices[2];
	triangle[2] = _indices[3];
	_faces.push_back(new Triangle3(triangle));
	_faces.back()->elements().push_back(this);

	triangle[0] = _indices[2];
	triangle[1] = _indices[0];
	triangle[2] = _indices[3];
	_faces.push_back(new Triangle3(triangle));
	_faces.back()->elements().push_back(this);
}

Tetrahedron4::Tetrahedron4(const eslocal *indices, eslocal n, const eslocal *params)
{
	switch (n) {
	case 8:
		memcpy(_indices, indices, 3 * sizeof(eslocal));
		_indices[3] = indices[4];
		break;
	case 4:
		memcpy(_indices, indices, 4 * sizeof(eslocal));
		break;
	default:
		ESINFO(ERROR) << "It is not possible to create Tetrahedron4 from " << n << " elements.";
	}

	memcpy(_params, params, PARAMS_SIZE * sizeof(eslocal));
}

Tetrahedron4::Tetrahedron4(std::ifstream &is)
{
	is.read(reinterpret_cast<char *>(_indices), sizeof(eslocal) * nodes());
	is.read(reinterpret_cast<char *>(_params), sizeof(eslocal) * PARAMS_SIZE);
}

