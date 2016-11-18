
#ifndef SRC_MESH_ELEMENTS_PLANE_TRIANGLE3_H_
#define SRC_MESH_ELEMENTS_PLANE_TRIANGLE3_H_

#include "planeelement.h"

#define Triangle3NodesCount 3
#define Triangle3EdgeCount 3
#define Triangle3GPCount 1
#define Triangle3CommonNodes 2
#define Triangle3VTKCode 5

namespace espreso {

class Triangle3: public PlaneElement
{

public:
	static bool match(eslocal *indices, eslocal n);
	static void setDOFs(
			const std::vector<Property> element,
			const std::vector<Property> face,
			const std::vector<Property> edge,
			const std::vector<Property> point,
			const std::vector<Property> midPoint)
	{
		_DOFElement = element;
		_DOFFace = face;
		_DOFEdge = edge;
		_DOFPoint = point;
		_DOFMidPoint = midPoint;
	}

	Triangle3(const eslocal *indices);
	Triangle3(const eslocal *indices, const eslocal *params);
	Triangle3(std::ifstream &is);
	Element* copy() const { return new Triangle3(*this); }

	eslocal nCommon() const { return Triangle3CommonNodes; }
	eslocal vtkCode() const { return Triangle3VTKCode; }

	size_t edges() const { return Triangle3EdgeCount; }
	size_t nodes() const { return Triangle3NodesCount; }
	size_t coarseNodes() const { return Triangle3NodesCount; }
	size_t gaussePoints() const { return Triangle3GPCount; }

	virtual Point edgeNormal(const Element *edge, const Coordinates &coordinates) const;

	const std::vector<DenseMatrix>& dN() const { return Triangle3::_dN; }
	const std::vector<DenseMatrix>& N() const { return Triangle3::_N; }
	const std::vector<double>& weighFactor() const { return Triangle3::_weighFactor; }

	const std::vector<Property>& elementDOFs() const { return Triangle3::_DOFElement; }
	const std::vector<Property>& faceDOFs() const { return Triangle3::_DOFFace; }
	const std::vector<Property>& edgeDOFs() const { return Triangle3::_DOFEdge; }
	const std::vector<Property>& pointDOFs() const { return Triangle3::_DOFPoint; }
	const std::vector<Property>& midPointDOFs() const { return Triangle3::_DOFMidPoint; }

protected:
	std::vector<eslocal> getNeighbours(size_t nodeIndex) const;
	eslocal* indices() { return _indices; }
	const eslocal* indices() const { return _indices; }

	void fillEdges();

private:
	eslocal _indices[Triangle3NodesCount];

	static std::vector<DenseMatrix> _dN;
	static std::vector<DenseMatrix> _N;
	static std::vector<double> _weighFactor;

	static std::vector<Property> _DOFElement;
	static std::vector<Property> _DOFFace;
	static std::vector<Property> _DOFEdge;
	static std::vector<Property> _DOFPoint;
	static std::vector<Property> _DOFMidPoint;
};

}



#endif /* SRC_MESH_ELEMENTS_PLANE_TRIANGLE3_H_ */
