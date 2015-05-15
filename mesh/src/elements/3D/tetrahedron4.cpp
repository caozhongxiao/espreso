
#include "tetrahedron4.h"

std::vector<std::vector<double> > Tetra4_dN()
{
	std::vector<std::vector<double> > dN(Tetrahedron4GPCount);

	// dN contains [dNr, dNs, dNt]
	int dN_length  = 12;

	for (unsigned int i = 0; i < Tetrahedron4GPCount; i++) {
		//  N = [ r, s, t,  1 - r - s - t ];
		dN[i].assign(dN_length, 0);

		// dNr = [ 1, 0, 0, -1 ];
		dN[i][0] =  1.0;
		dN[i][1] =  0.0;
		dN[i][2] =  0.0;
		dN[i][3] = -1.0;

		// dNs = [ 0, 1, 0, -1 ];
		dN[i][4] =  0.0;
		dN[i][5] =  1.0;
		dN[i][6] =  0.0;
		dN[i][7] = -1.0;

		// dNs = [ 0, 0, 1, -1 ];
		dN[i][8] =  0.0;
		dN[i][9] =  0.0;
		dN[i][10] =  1.0;
		dN[i][11] = -1.0;
	}

	return dN;
}

std::vector<std::vector<double> > Tetra4_N()
{
	std::vector<std::vector<double> > N(Tetrahedron4GPCount);

	std::vector<double> rv;
	std::vector<double> sv;
	std::vector<double> tv;

	if (Tetrahedron4GPCount == 4) {
		double _rv[] = {0.5854101966249685, 0.1381966011250105, 0.1381966011250105, 0.1381966011250105};
		double _sv[] = {0.1381966011250105, 0.1381966011250105, 0.1381966011250105, 0.5854101966249685};
		double _tv[] = {0.1381966011250105, 0.1381966011250105, 0.5854101966249685, 0.1381966011250105};
		rv.assign(_rv, _rv + Tetrahedron4GPCount);
		rv.assign(_sv, _sv + Tetrahedron4GPCount);
		rv.assign(_tv, _tv + Tetrahedron4GPCount);
	}
	else if (Tetrahedron4GPCount == 5) {
		double _rv[] = {0.2500000000000000, 0.5000000000000000, 0.1666666666666667, 0.1666666666666667, 0.1666666666666667};
		double _sv[] = {0.2500000000000000, 0.1666666666666667, 0.1666666666666667, 0.1666666666666667, 0.5000000000000000};
		double _tv[] = {0.2500000000000000, 0.1666666666666667, 0.1666666666666667, 0.5000000000000000, 0.1666666666666667};
		rv.assign(_rv, _rv + Tetrahedron4GPCount);
		rv.assign(_sv, _sv + Tetrahedron4GPCount);
		rv.assign(_tv, _tv + Tetrahedron4GPCount);
	}
	else if (Tetrahedron4GPCount == 11) {
		double _rv[] = {0.2500000000000000, 0.7857142857142857, 0.0714285714285714, 0.0714285714285714,
		0.0714285714285714, 0.1005964238332008, 0.3994035761667992, 0.3994035761667992,
		0.3994035761667992, 0.1005964238332008, 0.1005964238332008};
		double _sv[] = {0.2500000000000000, 0.0714285714285714, 0.0714285714285714, 0.0714285714285714,
		0.7857142857142857, 0.3994035761667992, 0.1005964238332008, 0.3994035761667992,
		0.1005964238332008, 0.3994035761667992, 0.1005964238332008};
		double _tv[] = {0.2500000000000000, 0.0714285714285714, 0.0714285714285714, 0.7857142857142857,
		0.0714285714285714, 0.3994035761667992, 0.3994035761667992, 0.1005964238332008,
		0.1005964238332008, 0.1005964238332008, 0.3994035761667992};
		rv.assign(_rv, _rv + Tetrahedron4GPCount);
		rv.assign(_sv, _sv + Tetrahedron4GPCount);
		rv.assign(_tv, _tv + Tetrahedron4GPCount);
	}
	else if (Tetrahedron4GPCount == 15) {
		double _rv[] = {0.2500000000000000, 0.0000000000000000, 0.3333333333333333, 0.3333333333333333,
		0.3333333333333333, 0.7272727272727273, 0.0909090909090909, 0.0909090909090909,
		0.0909090909090909, 0.4334498464263357, 0.0665501535736643, 0.0665501535736643,
		0.0665501535736643, 0.4334498464263357, 0.4334498464263357};
		double _sv[] = {0.2500000000000000, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
		0.0000000000000000, 0.0909090909090909, 0.0909090909090909, 0.0909090909090909,
		0.7272727272727273, 0.0665501535736643, 0.4334498464263357, 0.0665501535736643,
		0.4334498464263357, 0.0665501535736643, 0.4334498464263357};
		double _tv[] = {0.2500000000000000, 0.3333333333333333, 0.3333333333333333, 0.0000000000000000,
		0.3333333333333333, 0.0909090909090909, 0.0909090909090909, 0.7272727272727273,
		0.0909090909090909, 0.0665501535736643, 0.0665501535736643, 0.4334498464263357,
		0.4334498464263357, 0.4334498464263357, 0.0665501535736643};
		rv.assign(_rv, _rv + Tetrahedron4GPCount);
		rv.assign(_sv, _sv + Tetrahedron4GPCount);
		rv.assign(_tv, _tv + Tetrahedron4GPCount);
	}


	// N = [r s t  1-r-s-t];
	for (unsigned int i = 0; i < Tetrahedron4GPCount; i++) {
		double r = rv[i];
		double s = sv[i];
		double t = tv[i];

		N[i].resize(4);
		N[i][0] = r;
		N[i][1] = s;
		N[i][2] = t;
		N[i][3] = 1.0 - r - s - t;
	}

	return N;
}

std::vector<double> Tetra4_Weight()
{
	switch (Tetrahedron4GPCount) {
	case 4: {
		return std::vector<double> (4, 1 / 24.0);
	}
	case 5: {
		std::vector<double> w(5, 3 / 40.0);
		w[0] = - 2 / 15.0;
		return w;
	}
	case 11: {
		std::vector<double> w(11);
		w[0] = -0.013155555555556;
		w[1] = w[2] = w[3] = w[4] = 0.007622222222222;
		w[5] = w[6] = w[7] = w[8] = w[5] = w[6] = w[7] = 0.024888888888889;
		return w;
	}
	case 15: {
		std::vector<double> w(15);
		w[0] = 0.030283678097089;
		w[1] = w[2] = w[3] = w[4] = 0.006026785714286;
		w[5] = w[6] = w[7] = w[8] = w[9] = w[10] = w[5] = w[6] = w[7] = w[8] = w[9] = 0.011645249086029;
		return w;
	}
	default:
		std::cerr << "Unknown number of Tatrahedron4 GP count\n";
		exit(EXIT_FAILURE);
	}
}

std::vector<std::vector<double> > Tetrahedron4::_dN = Tetra4_dN();
std::vector<std::vector<double> > Tetrahedron4::_N = Tetra4_N();
std::vector<double> Tetrahedron4::_weighFactor = Tetra4_Weight();

bool Tetrahedron4::match(idx_t *indices, idx_t n) {

#ifndef D3
	// Tetrahedron4 is 3D element
	return false;
#endif

	if (n != 8) {
		return false;
	}

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

	idx_t various[4] = { 0, 1, 2, 4 };
	for (idx_t i = 0; i < 3; i++) {
		for (idx_t j = i + 1; j < 4; j++) {
			if (Element::match(indices, various[i], various[j])) {
				return false;
			}
		}
	}

	return true;
}

std::vector<idx_t> Tetrahedron4::getNeighbours(size_t nodeIndex) const
{
	std::vector<idx_t> result;
	result.reserve(3);
	for (size_t i = 0; i < Tetrahedron4NodesCount; i++) {
		if (i != nodeIndex) {
			result.push_back(_indices[i]);
		}
	}
	return result;
}

std::vector<idx_t> Tetrahedron4::getFace(size_t face) const
{
	std::vector<idx_t> result(3);
	result[0] = (face < 3) ? _indices[0] : _indices[1];
	result[1] = (face < 2) ? _indices[1] : _indices[2];
	result[2] = (face < 1) ? _indices[2] : _indices[3];
	return result;
}

Tetrahedron4::Tetrahedron4(idx_t *indices)
{
	memcpy(_indices, indices, Tetrahedron4NodesCount * sizeof(idx_t));
}


