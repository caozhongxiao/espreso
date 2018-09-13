
#include "prisma6.h"

using namespace espreso;

Prisma6Generator::Prisma6Generator()
{
	subelements = 2;
	enodes = 6;
	code = Element::CODE::PRISMA6;
}

void Prisma6Generator::pushElements(std::vector<eslocal> &elements, const std::vector<eslocal> &indices) const
{
	elements.push_back(indices[0]);
	elements.push_back(indices[1]);
	elements.push_back(indices[3]);
	elements.push_back(indices[4]);
	elements.push_back(indices[5]);
	elements.push_back(indices[7]);

	elements.push_back(indices[0]);
	elements.push_back(indices[3]);
	elements.push_back(indices[2]);
	elements.push_back(indices[4]);
	elements.push_back(indices[7]);
	elements.push_back(indices[6]);
}

void Prisma6Generator::pushNodes(std::vector<eslocal> &nodes, const std::vector<eslocal> &indices, CubeFace face) const
{
	switch (face) {
	case CubeFace::X_0:
	case CubeFace::X_1:
	case CubeFace::Y_0:
	case CubeFace::Y_1:
		pushSquareNodes(nodes, indices, face);
		break;
	case CubeFace::Z_0:
	case CubeFace::Z_1:
		pushTriangleNodes(nodes, indices, face);
		break;
	}
}

void Prisma6Generator::pushFace(std::vector<eslocal> &elements, std::vector<eslocal> &esize, std::vector<int> &etype, const std::vector<eslocal> &indices, CubeFace face) const
{
	pushNodes(elements, indices, face);
	switch (face) {
	case CubeFace::X_0:
	case CubeFace::X_1:
	case CubeFace::Y_0:
	case CubeFace::Y_1:
		esize.push_back(4);
		etype.push_back((int)Element::CODE::SQUARE4);
		break;
	case CubeFace::Z_0:
	case CubeFace::Z_1:
		esize.push_back(3);
		esize.push_back(3);
		etype.push_back((int)Element::CODE::TRIANGLE3);
		etype.push_back((int)Element::CODE::TRIANGLE3);
		break;
	}
}
