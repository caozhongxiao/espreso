
#ifndef SRC_MESH_ELEMENTS_VOLUME_PRISMA15_H_
#define SRC_MESH_ELEMENTS_VOLUME_PRISMA15_H_

#include "volumeelement.h"
#include "../line/line3.h"
#include "../plane/square8.h"
#include "../plane/triangle6.h"
#include "prisma6.h"

#define Prisma15NodesCount 15
#define Prisma15EdgeCount 9
#define Prisma15FacesCount 5
#define Prisma15GPCount 9
#define Prisma15CommonNodes 3
#define Prisma15VTKCode 26

namespace espreso {

class Prisma15: public VolumeElement
{

public:
	static bool match(const eslocal *indices, eslocal n);
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

	Prisma15(const eslocal *indices, eslocal n, const eslocal *params);
	Prisma15(std::ifstream &is);
	OldElement* copy() const { return new Prisma15(*this); }

	eslocal nCommon() const { return Prisma15CommonNodes; }
	eslocal vtkCode() const { return Prisma15VTKCode; }

	size_t faces() const { return Prisma15FacesCount; }
	size_t edges() const { return Prisma15EdgeCount; }
	size_t nodes() const { return Prisma15NodesCount; }
	size_t coarseNodes() const { return Prisma6NodesCount; }
	size_t gaussePoints() const { return Prisma15GPCount; }

	OldElement* addFace(const std::vector<eslocal> &nodes);

	const std::vector<eslocal>& faceNodes(size_t index) const { return Prisma15::_facesNodes[index]; }
	const std::vector<eslocal>& edgeNodes(size_t index) const { return Prisma15::_edgesNodes[index]; }

	const std::vector<DenseMatrix>& facedN(size_t index, ElementPointType type = ElementPointType::GAUSSE_POINT) const { return index < 3 ? Square8::_dN : Triangle6::_dN; }
	const std::vector<DenseMatrix>& faceN(size_t index, ElementPointType type = ElementPointType::GAUSSE_POINT) const { return index < 3 ? Square8::_N : Triangle6::_N; }
	const std::vector<DenseMatrix>& edgedN(size_t index, ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Line3::_dN; }
	const std::vector<DenseMatrix>& edgeN(size_t index, ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Line3::_N; }

	const std::vector<DenseMatrix>& dN(ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Prisma15::_dN; }
	const std::vector<DenseMatrix>& N(ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Prisma15::_N; }
	const std::vector<double>& weighFactor(ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Prisma15::_weighFactor; }

	const std::vector<Property>& elementDOFs() const { return Prisma15::_DOFElement; }
	const std::vector<Property>& faceDOFs() const { return Prisma15::_DOFFace; }
	const std::vector<Property>& edgeDOFs() const { return Prisma15::_DOFEdge; }
	const std::vector<Property>& pointDOFs() const { return Prisma15::_DOFPoint; }
	const std::vector<Property>& midPointDOFs() const { return Prisma15::_DOFMidPoint; }

	static std::vector<DenseMatrix> _dN;
	static std::vector<DenseMatrix> _N;
	static std::vector<double> _weighFactor;

protected:
	std::vector<eslocal> getNeighbours(size_t nodeIndex) const;
	eslocal* indices() { return _indices; }
	const eslocal* indices() const { return _indices; }

	size_t fillFaces();
	size_t fillEdges();

private:
	eslocal _indices[Prisma15NodesCount];

	static std::vector<Property> _DOFElement;
	static std::vector<Property> _DOFFace;
	static std::vector<Property> _DOFEdge;
	static std::vector<Property> _DOFPoint;
	static std::vector<Property> _DOFMidPoint;

	static std::vector<std::vector<eslocal> > _facesNodes;
	static std::vector<std::vector<eslocal> > _edgesNodes;
};

}


#endif /* SRC_MESH_ELEMENTS_VOLUME_PRISMA15_H_ */
