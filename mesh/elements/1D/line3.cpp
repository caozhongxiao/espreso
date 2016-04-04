#include "line3.h"

using namespace espreso;

// TODO: Implement base functions
std::vector<DenseMatrix> Line3::_dN;
std::vector<DenseMatrix> Line3::_N;
std::vector<double> Line3::_weighFactor;

bool Line3::match(const eslocal *indices, eslocal n)
{
	if (n != 2) {
		return false;
	}

	if (Element::match(indices, 0, 1)) {
		return false;
	}

	return true;
}

std::vector<eslocal> Line3::getNeighbours(size_t nodeIndex) const
{
	std::vector<eslocal> result(1);

	result[0] = (nodeIndex == 0) ? _indices[1] : _indices[0];

	return result;
}

std::vector<eslocal> Line3::getFace(size_t face) const
{
	return std::vector<eslocal> (0);
}

Element* Line3::getFullFace(size_t face) const
{
	ESINFO(ERROR) << "get FACE is not possible on Line3";
	return NULL;
}

Element* Line3::getCoarseFace(size_t face) const
{
	ESINFO(ERROR) << "get FACE is not possible on Line3";
	return NULL;
}

Line3::Line3(const eslocal *indices, const eslocal *params): Element(params)
{
	memcpy(_indices, indices, Line3NodesCount * sizeof(eslocal));
}



