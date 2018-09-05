
#include "linearplane.h"

using namespace espreso;

LinearPlaneGenerator::LinearPlaneGenerator()
{
	subnodes[0] = 2;
	subnodes[1] = 2;
	subnodes[2] = 1;
}

void LinearPlaneGenerator::pushNodes(std::vector<eslocal> &nodes, const std::vector<eslocal> &indices, CubeEdge edge) const
{
	switch (edge) {
	case CubeEdge::X_0_Z_0:
		nodes.push_back(indices[2]);
		nodes.push_back(indices[0]);
		break;
	case CubeEdge::X_1_Z_0:
		nodes.push_back(indices[1]);
		nodes.push_back(indices[3]);
		break;
	case CubeEdge::Y_0_Z_0:
		nodes.push_back(indices[0]);
		nodes.push_back(indices[1]);
		break;
	case CubeEdge::Y_1_Z_0:
		nodes.push_back(indices[3]);
		nodes.push_back(indices[2]);
		break;
	default:
		return;
	}
}

void LinearPlaneGenerator::pushEdge(std::vector<eslocal> &elements, std::vector<eslocal> &esize, std::vector<eslocal> &etype, const std::vector<eslocal> &indices, CubeEdge edge) const
{
	pushNodes(elements, indices, edge);
	esize.push_back(2);
	etype.push_back((eslocal)Element::CODE::LINE2);
}

void LinearPlaneGenerator::pushNodes(std::vector<eslocal> &nodes, const std::vector<eslocal> &indices, CubeFace face) const
{
	return;
}

void LinearPlaneGenerator::pushFace(std::vector<eslocal> &elements, std::vector<eslocal> &esize, std::vector<eslocal> &etype, const std::vector<eslocal> &indices, CubeFace face) const
{
	return;
}


