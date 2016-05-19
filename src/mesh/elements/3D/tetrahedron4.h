
#ifndef TETRAHEDRON4_H_
#define TETRAHEDRON4_H_

#include "../2D/triangle3.h"
#include "../element.h"

#define Tetrahedron4NodesCount 4
#define Tetrahedron4FacesCount 4
#define Tetrahedron4GPCount 4
#define Tetrahedron4VTKCode 10

namespace espreso {

class Tetrahedron4: public Element
{
	friend class Tetrahedron10;

public:
	static bool match(const eslocal *indices, eslocal n);

	Tetrahedron4(const eslocal *indices, eslocal n, const eslocal *params);
	Tetrahedron4(std::ifstream &is);

	Element* copy() const
	{
		return new Tetrahedron4(*this);
	}

	eslocal vtkCode() const
	{
		return Tetrahedron4VTKCode;
	}

	const eslocal* indices() const
	{
		return _indices;
	}

	size_t size() const
	{
		return Tetrahedron4NodesCount;
	}

	size_t coarseSize() const
	{
		return Tetrahedron4NodesCount;
	}

	size_t gpSize() const
	{
		return Tetrahedron4GPCount;
	}

	size_t faces() const
	{
		return Tetrahedron4FacesCount;
	}

	const std::vector<DenseMatrix>& dN() const
	{
		return Tetrahedron4::_dN;
	}

	const std::vector<DenseMatrix>& N() const
	{
		return Tetrahedron4::_N;
	}

	const std::vector<double>& weighFactor() const
	{
		return Tetrahedron4::_weighFactor;
	}

	eslocal nCommon() const
	{
		return 3;
	}

	Element* getFullFace(size_t face) const
	{
		return getF(_indices, _params, face);
	}

	Element* getCoarseFace(size_t face) const
	{
		return getFullFace(face);
	}

	std::vector<eslocal> getNeighbours(size_t nodeIndex) const;
	std::vector<eslocal> getFace(size_t face) const;

protected:
	static Element* getF(const eslocal *indices, const eslocal *params, size_t face);

	eslocal* indices()
	{
		return _indices;
	}

private:
	eslocal _indices[Tetrahedron4NodesCount];

	static std::vector<DenseMatrix> _dN;
	static std::vector<DenseMatrix> _N;
	static std::vector<double> _weighFactor;
};


}

#endif /* TETRAHEDRON4_H_ */