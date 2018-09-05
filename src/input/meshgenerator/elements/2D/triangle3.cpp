
#include "triangle3.h"

using namespace espreso;

Triangle3Generator::Triangle3Generator()
{
	subelements = 2;
	enodes = 3;
	code = Element::CODE::TRIANGLE3;
}

void Triangle3Generator::pushElements(std::vector<eslocal> &elements, const std::vector<eslocal> &indices) const
{
	elements.push_back(indices[0]);
	elements.push_back(indices[1]);
	elements.push_back(indices[3]);

	elements.push_back(indices[0]);
	elements.push_back(indices[3]);
	elements.push_back(indices[2]);
}



