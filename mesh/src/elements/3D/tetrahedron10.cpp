#include "tetrahedron10.h"

std::vector<std::vector<double> > Tetra10_dN()
{
	std::vector<std::vector<double> > dN(Tetrahedron10GPCount);

	int dN_length  = 30;

	std::vector<double> rv;
	std::vector<double> sv;
	std::vector<double> tv;



	if (Tetrahedron10GPCount == 4) {
		double _rv[] = {0.5854101966249685, 0.1381966011250105, 0.1381966011250105, 0.1381966011250105};
		double _sv[] = {0.1381966011250105, 0.1381966011250105, 0.1381966011250105, 0.5854101966249685};
		double _tv[] = {0.1381966011250105, 0.1381966011250105, 0.5854101966249685, 0.1381966011250105};
		rv.assign(_rv, _rv + Tetrahedron10GPCount);
		sv.assign(_sv, _sv + Tetrahedron10GPCount);
		tv.assign(_tv, _tv + Tetrahedron10GPCount);
	}
	else if (Tetrahedron10GPCount == 5) {
		double _rv[] = {0.2500000000000000, 0.5000000000000000, 0.1666666666666667, 0.1666666666666667, 0.1666666666666667};
		double _sv[] = {0.2500000000000000, 0.1666666666666667, 0.1666666666666667, 0.1666666666666667, 0.5000000000000000};
		double _tv[] = {0.2500000000000000, 0.1666666666666667, 0.1666666666666667, 0.5000000000000000, 0.1666666666666667};
		rv.assign(_rv, _rv + Tetrahedron10GPCount);
		sv.assign(_sv, _sv + Tetrahedron10GPCount);
		tv.assign(_tv, _tv + Tetrahedron10GPCount);
	}
	else if (Tetrahedron10GPCount == 11) {
		double _rv[] = {0.2500000000000000, 0.7857142857142857, 0.0714285714285714, 0.0714285714285714,
		0.0714285714285714, 0.1005964238332008, 0.3994035761667992, 0.3994035761667992,
		0.3994035761667992, 0.1005964238332008, 0.1005964238332008};
		double _sv[] = {0.2500000000000000, 0.0714285714285714, 0.0714285714285714, 0.0714285714285714,
		0.7857142857142857, 0.3994035761667992, 0.1005964238332008, 0.3994035761667992,
		0.1005964238332008, 0.3994035761667992, 0.1005964238332008};
		double _tv[] = {0.2500000000000000, 0.0714285714285714, 0.0714285714285714, 0.7857142857142857,
		0.0714285714285714, 0.3994035761667992, 0.3994035761667992, 0.1005964238332008,
		0.1005964238332008, 0.1005964238332008, 0.3994035761667992};
		rv.assign(_rv, _rv + Tetrahedron10GPCount);
		sv.assign(_sv, _sv + Tetrahedron10GPCount);
		tv.assign(_tv, _tv + Tetrahedron10GPCount);
	}
	else if (Tetrahedron10GPCount == 15) {
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
		rv.assign(_rv, _rv + Tetrahedron10GPCount);
		sv.assign(_sv, _sv + Tetrahedron10GPCount);
		tv.assign(_tv, _tv + Tetrahedron10GPCount);
	}

	for (unsigned int i = 0; i < Tetrahedron10GPCount; i++) {
		double r = rv[i];
		double s = sv[i];
		double t = tv[i];

		dN[i].assign(dN_length, 0);

//	N[0] = r*(2.0*r - 1.0) ;
//	N[1] = s*(2.0*s - 1.0) ;
//	N[2] = t*(2.0*t - 1.0) ;
//	N[3] = 2.0*r**2 + 4.0*r*s + 4.0*r*t - 3.0*r + 2.0*s**2 + 4.0*s*t - 3.0*s + 2.0*t**2 - 3.0*t + 1.0 ;
//	N[4] = 4.0*r*s ;
//	N[5] = 4.0*s*t ;
//	N[6] = 4.0*r*t ;
//	N[7] = r*(-4.0*r - 4.0*s - 4.0*t + 4.0) ;
//	N[8] = s*(-4.0*r - 4.0*s - 4.0*t + 4.0) ;
//	N[9] = t*(-4.0*r - 4.0*s - 4.0*t + 4.0) ;


    dN[i][ 0] = 4.0*r - 1.0;
    dN[i][ 1] = 0;
    dN[i][ 2] = 0;
    dN[i][ 3] = 4.0*r + 4.0*s + 4.0*t - 3.0;
    dN[i][ 4] = 4.0*s;
    dN[i][ 5] = 0;
    dN[i][ 6] = 4.0*t;
    dN[i][ 7] = -8.0*r - 4.0*s - 4.0*t + 4.0;
    dN[i][ 8] = -4.0*s;
    dN[i][ 9] = -4.0*t;

    dN[i][10] = 0 ;
    dN[i][11] = 4.0*s - 1.0 ;
    dN[i][12] = 0 ;
    dN[i][13] = 4.0*r + 4.0*s + 4.0*t - 3.0 ;
    dN[i][14] = 4.0*r ;
    dN[i][15] = 4.0*t ;
    dN[i][16] = 0 ;
    dN[i][17] = -4.0*r ;
    dN[i][18] = -4.0*r - 8.0*s - 4.0*t + 4.0 ;
    dN[i][19] = -4.0*t ;

    dN[i][20] = 0 ;
    dN[i][21] = 0 ;
    dN[i][22] = 4.0*t - 1.0 ;
    dN[i][23] = 4.0*r + 4.0*s + 4.0*t - 3.0 ;
    dN[i][24] = 0 ;
    dN[i][25] = 4.0*s ;
    dN[i][26] = 4.0*r ;
    dN[i][27] = -4.0*r ;
    dN[i][28] = -4.0*s ;
    dN[i][29] = -4.0*r - 4.0*s - 8.0*t + 4.0 ;

	}

	return dN;
}

std::vector<std::vector<double> > Tetra10_N() {
	std::vector<std::vector<double> > N(Tetrahedron10GPCount);

	std::vector<double> rv;
	std::vector<double> sv;
	std::vector<double> tv;


	if (Tetrahedron10GPCount == 4) {
		double _rv[] = {0.5854101966249685, 0.1381966011250105, 0.1381966011250105, 0.1381966011250105};
		double _sv[] = {0.1381966011250105, 0.1381966011250105, 0.1381966011250105, 0.5854101966249685};
		double _tv[] = {0.1381966011250105, 0.1381966011250105, 0.5854101966249685, 0.1381966011250105};
		rv.assign(_rv, _rv + Tetrahedron10GPCount);
		sv.assign(_sv, _sv + Tetrahedron10GPCount);
		tv.assign(_tv, _tv + Tetrahedron10GPCount);
	}
	else if (Tetrahedron10GPCount == 5) {
		double _rv[] = {0.2500000000000000, 0.5000000000000000, 0.1666666666666667, 0.1666666666666667, 0.1666666666666667};
		double _sv[] = {0.2500000000000000, 0.1666666666666667, 0.1666666666666667, 0.1666666666666667, 0.5000000000000000};
		double _tv[] = {0.2500000000000000, 0.1666666666666667, 0.1666666666666667, 0.5000000000000000, 0.1666666666666667};
		rv.assign(_rv, _rv + Tetrahedron10GPCount);
		sv.assign(_sv, _sv + Tetrahedron10GPCount);
		tv.assign(_tv, _tv + Tetrahedron10GPCount);
	}
	else if (Tetrahedron10GPCount == 11) {
		double _rv[] = {0.2500000000000000, 0.7857142857142857, 0.0714285714285714, 0.0714285714285714,
		0.0714285714285714, 0.1005964238332008, 0.3994035761667992, 0.3994035761667992,
		0.3994035761667992, 0.1005964238332008, 0.1005964238332008};
		double _sv[] = {0.2500000000000000, 0.0714285714285714, 0.0714285714285714, 0.0714285714285714,
		0.7857142857142857, 0.3994035761667992, 0.1005964238332008, 0.3994035761667992,
		0.1005964238332008, 0.3994035761667992, 0.1005964238332008};
		double _tv[] = {0.2500000000000000, 0.0714285714285714, 0.0714285714285714, 0.7857142857142857,
		0.0714285714285714, 0.3994035761667992, 0.3994035761667992, 0.1005964238332008,
		0.1005964238332008, 0.1005964238332008, 0.3994035761667992};
		rv.assign(_rv, _rv + Tetrahedron10GPCount);
		sv.assign(_sv, _sv + Tetrahedron10GPCount);
		tv.assign(_tv, _tv + Tetrahedron10GPCount);
	}
	else if (Tetrahedron10GPCount == 15) {
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
		rv.assign(_rv, _rv + Tetrahedron10GPCount);
		sv.assign(_sv, _sv + Tetrahedron10GPCount);
		tv.assign(_tv, _tv + Tetrahedron10GPCount);
	}

	for (unsigned int i = 0; i < Tetrahedron10GPCount; i++) {
		double r = rv[i];
		double s = sv[i];
		double t = tv[i];

		N[i].resize(10);


    N[i][0] = r*(2.0*r - 1.0) ;
    N[i][1] = s*(2.0*s - 1.0) ;
    N[i][2] = t*(2.0*t - 1.0) ;
    N[i][3] = 2.0*r*r + 4.0*r*s + 4.0*r*t - 3.0*r + 2.0*s*s + 4.0*s*t - 3.0*s + 2.0*t*t - 3.0*t + 1.0 ;
    N[i][4] = 4.0*r*s ;
    N[i][5] = 4.0*s*t ;
    N[i][6] = 4.0*r*t ;
    N[i][7] = r*(-4.0*r - 4.0*s - 4.0*t + 4.0) ;
    N[i][8] = s*(-4.0*r - 4.0*s - 4.0*t + 4.0) ;
    N[i][9] = t*(-4.0*r - 4.0*s - 4.0*t + 4.0) ;

	}

	return N;
}



std::vector<double> Tetra10_Weight()
{
	switch (Tetrahedron10GPCount) {
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
		std::cerr << "Unknown number of Tatrahedron10 GP count\n";
		exit(EXIT_FAILURE);
	}
}



std::vector<std::vector<double> > Tetrahedron10::_dN = Tetra10_dN();
std::vector<std::vector<double> > Tetrahedron10::_N = Tetra10_N();
std::vector<double> Tetrahedron10::_weighFactor = Tetra10_Weight();

bool Tetrahedron10::match(idx_t *indices, idx_t n) {

#ifndef D3
	// Tetrahedron10 is 3D element
	return false;
#endif

	if (n != 20) {
		return false;
	}

	if (!Element::match(indices, 2, 3)) {
		return false;
	}
	if (!Element::match(indices, 2, 10)) {
		return false;
	}
	if (!Element::match(indices, 19, 20)) {
		return false;
	}
	if (!Element::match(indices, 4, 5)) {
		return false;
	}
	if (!Element::match(indices, 4, 6)) {
		return false;
	}
	if (!Element::match(indices, 4, 7)) {
		return false;
	}
	if (!Element::match(indices, 4, 12)) {
		return false;
	}
	if (!Element::match(indices, 4, 13)) {
		return false;
	}
	if (!Element::match(indices, 4, 14)) {
		return false;
	}
	if (!Element::match(indices, 4, 15)) {
		return false;
	}

	idx_t various[10] = { 0, 1, 2, 4, 8, 9 ,11, 16, 17, 18 };
	for (idx_t i = 0; i < 9; i++) {
		for (idx_t j = i + 1; j < 10; j++) {
			if (Element::match(indices, various[i], various[j])) {
				return false;
			}
		}
	}

	return true;
}

std::vector<idx_t> Tetrahedron10::getNeighbours(size_t nodeIndex) const
{
	std::vector<idx_t> result;
	if (nodeIndex > 3) {
		result.resize(2);
	} else {
		result.resize(3);
	}

	switch (nodeIndex) {
	case 0: {
		result[0] = _indices[4];
		result[1] = _indices[6];
		result[2] = _indices[7];
		return result;
	}
	case 1: {
		result[0] = _indices[4];
		result[1] = _indices[5];
		result[2] = _indices[8];
		return result;
	}
	case 2: {
		result[0] = _indices[5];
		result[1] = _indices[6];
		result[2] = _indices[9];
		return result;
	}
	case 3: {
		result[0] = _indices[7];
		result[1] = _indices[8];
		result[2] = _indices[9];
		return result;
	}
	case 4: {
		result[0] = _indices[0];
		result[1] = _indices[1];
		return result;
	}
	case 5: {
		result[0] = _indices[1];
		result[1] = _indices[2];
		return result;
	}
	case 6: {
		result[0] = _indices[0];
		result[1] = _indices[2];
		return result;
	}
	case 7: {
		result[0] = _indices[0];
		result[1] = _indices[3];
		return result;
	}
	case 8: {
		result[0] = _indices[1];
		result[1] = _indices[3];
		return result;
	}
	case 9: {
		result[0] = _indices[2];
		result[1] = _indices[3];
		return result;
	}
	}
	return result;
}

std::vector<idx_t> Tetrahedron10::getFace(size_t face) const
{
	std::vector<idx_t> result(3);
	result[0] = (face < 3) ? _indices[0] : _indices[1];
	result[1] = (face < 2) ? _indices[1] : _indices[2];
	result[2] = (face < 1) ? _indices[2] : _indices[3];
	return result;
}

Tetrahedron10::Tetrahedron10(idx_t *indices)
{
	_indices[0] = indices[0];
	_indices[1] = indices[1];
	_indices[2] = indices[2];
	_indices[3] = indices[4];
	_indices[4] = indices[8];
	_indices[5] = indices[9];
	_indices[6] = indices[11];
	_indices[7] = indices[16];
	_indices[8] = indices[17];
	_indices[9] = indices[18];
}



