
#include "tetrahedron4.h"

using namespace espreso::input;

size_t Tetrahedron4::subelements = 6;
size_t Tetrahedron4::subnodes[] = { 0, 0, 0 };

void Tetrahedron4::addElements(std::vector<Element*> &elements, const eslocal indices[], const eslocal params[])
{
	eslocal tetra[4];
	tetra[0] = indices[0];
	tetra[1] = indices[3];
	tetra[2] = indices[2];
	tetra[3] = indices[4];
	elements.push_back(new espreso::Tetrahedron4(tetra, 4, params));

	tetra[0] = indices[3];
	tetra[1] = indices[2];
	tetra[2] = indices[4];
	tetra[3] = indices[6];
	elements.push_back(new espreso::Tetrahedron4(tetra, 4, params));

	tetra[0] = indices[7];
	tetra[1] = indices[3];
	tetra[2] = indices[4];
	tetra[3] = indices[6];
	elements.push_back(new espreso::Tetrahedron4(tetra, 4, params));

	tetra[0] = indices[3];
	tetra[1] = indices[5];
	tetra[2] = indices[7];
	tetra[3] = indices[4];
	elements.push_back(new espreso::Tetrahedron4(tetra, 4, params));

	tetra[0] = indices[1];
	tetra[1] = indices[5];
	tetra[2] = indices[3];
	tetra[3] = indices[4];
	elements.push_back(new espreso::Tetrahedron4(tetra, 4, params));

	tetra[0] = indices[0];
	tetra[1] = indices[4];
	tetra[2] = indices[1];
	tetra[3] = indices[3];
	elements.push_back(new espreso::Tetrahedron4(tetra, 4, params));
}


