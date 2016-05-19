
#include "tetrahedron10.h"

using namespace espreso::input;

size_t Tetrahedron10::subelements = 6;
size_t Tetrahedron10::subnodes[] = { 1, 1, 1 };

void Tetrahedron10::addElements(std::vector<Element*> &elements, const eslocal indices[], const eslocal params[])
{
	eslocal size = 10;
	eslocal tetra[size];
	tetra[0] = indices[2];
	tetra[1] = indices[6];
	tetra[2] = indices[0];
	tetra[3] = indices[20];

	tetra[4] = indices[4];
	tetra[5] = indices[3];
	tetra[6] = indices[1];
	tetra[7] = indices[11];
	tetra[8] = indices[13];
	tetra[9] = indices[10];
	elements.push_back(new espreso::Tetrahedron10(tetra, size, params));

	tetra[0] = indices[6];
	tetra[1] = indices[0];
	tetra[2] = indices[20];
	tetra[3] = indices[18];

	tetra[4] = indices[3];
	tetra[5] = indices[10];
	tetra[6] = indices[13];
	tetra[7] = indices[12];
	tetra[8] = indices[9];
	tetra[9] = indices[19];
	elements.push_back(new espreso::Tetrahedron10(tetra, size, params));

	tetra[0] = indices[24];
	tetra[1] = indices[6];
	tetra[2] = indices[20];
	tetra[3] = indices[18];

	tetra[4] = indices[15];
	tetra[5] = indices[13];
	tetra[6] = indices[22];
	tetra[7] = indices[21];
	tetra[8] = indices[12];
	tetra[9] = indices[19];
	elements.push_back(new espreso::Tetrahedron10(tetra, size, params));

	tetra[0] = indices[6];
	tetra[1] = indices[26];
	tetra[2] = indices[24];
	tetra[3] = indices[20];

	tetra[4] = indices[16];
	tetra[5] = indices[25];
	tetra[6] = indices[15];
	tetra[7] = indices[13];
	tetra[8] = indices[23];
	tetra[9] = indices[22];
	elements.push_back(new espreso::Tetrahedron10(tetra, size, params));

	tetra[0] = indices[8];
	tetra[1] = indices[26];
	tetra[2] = indices[6];
	tetra[3] = indices[20];

	tetra[4] = indices[17];
	tetra[5] = indices[16];
	tetra[6] = indices[7];
	tetra[7] = indices[14];
	tetra[8] = indices[23];
	tetra[9] = indices[13];
	elements.push_back(new espreso::Tetrahedron10(tetra, size, params));

	tetra[0] = indices[2];
	tetra[1] = indices[20];
	tetra[2] = indices[8];
	tetra[3] = indices[6];

	tetra[4] = indices[11];
	tetra[5] = indices[14];
	tetra[6] = indices[5];
	tetra[7] = indices[4];
	tetra[8] = indices[13];
	tetra[9] = indices[7];
	elements.push_back(new espreso::Tetrahedron10(tetra, size, params));
}

