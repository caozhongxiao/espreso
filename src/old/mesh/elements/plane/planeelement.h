
#ifndef SRC_MESH_ELEMENTS_PLANE_PLANEELEMENT_H_
#define SRC_MESH_ELEMENTS_PLANE_PLANEELEMENT_H_

#include "../element.h"

namespace espreso {

class PlaneElement: public OldElement
{

public:
	Type type() const { return Type::PLANE; }

	virtual std::vector<std::vector<eslocal> > triangularize() const=0;

	eslocal param(Params param) const { return _params[param]; }
	void setParam(Params param, eslocal value) { _params[param] = value; }
	size_t params() const { return _params.size(); }

	virtual size_t filledFaces() const { return 0; }
	virtual size_t filledEdges() const { return _edges.size(); }

	size_t faces() const { return 0; }

	virtual OldElement* face(size_t index) const
	{
		ESINFO(GLOBAL_ERROR) << "Plane element has no face";
		return NULL;
	}
	virtual OldElement* edge(size_t index) const { return _edges[index]; }

	void addFace(OldElement* face) { ESINFO(GLOBAL_ERROR) << "Plane element has no face"; }
	virtual void addEdge(OldElement* edge)
	{
		_edges.push_back(edge);
	}

	OldElement* addFace(const std::vector<eslocal> &nodes) { return NULL; }

	const std::vector<eslocal>& faceNodes(size_t index) const
	{
		static std::vector<eslocal> _facesNodes;
		return _facesNodes;
	}

	const std::vector<DenseMatrix>& facedN(size_t index, ElementPointType type = ElementPointType::GAUSSE_POINT) const { ESINFO(ERROR) << "Plane element has no base functions for face."; exit(1); }
	const std::vector<DenseMatrix>& faceN(size_t index, ElementPointType type = ElementPointType::GAUSSE_POINT) const { ESINFO(ERROR) << "Plane element has no base functions for face."; exit(1); }

protected:
	void setFace(size_t index, OldElement* face) { ESINFO(GLOBAL_ERROR) << "Plane element has no face"; }
	void setEdge(size_t index, OldElement* edge) { _edges[index] = edge; }

	size_t fillFaces() { ESINFO(GLOBAL_ERROR) << "Call fill faces on plane element."; return 0; }

	std::vector<eslocal> _params;
	std::vector<OldElement*> _edges;
};

}

#endif /* SRC_MESH_ELEMENTS_PLANE_PLANEELEMENT_H_ */
