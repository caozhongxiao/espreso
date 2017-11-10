
#ifndef SRC_MESH_ELEMENTS_VOLUME_VOLUMEELEMENT_H_
#define SRC_MESH_ELEMENTS_VOLUME_VOLUMEELEMENT_H_

#include "../element.h"

namespace espreso {

class VolumeElement: public OldElement
{

	friend class OldMesh;

public:
	Type type() const { return Type::VOLUME; }

	virtual eslocal param(Params param) const { return _params[param]; };
	virtual void setParam(Params param, eslocal value) { _params[param] = value; }
	virtual size_t params() const { return PARAMS_SIZE; }

	virtual size_t filledFaces() const { return _faces.size(); }
	virtual size_t filledEdges() const { return _edges.size(); }

	virtual OldElement* face(size_t index) const
	{
		if (index >= _faces.size()) {
			ESINFO(ERROR) << "Element does not have face " << index;
		}
		return _faces[index];
	}
	virtual OldElement* edge(size_t index) const { return _edges[index]; }

	virtual void addFace(OldElement* face)
	{
		_faces.push_back(face);
	}
	virtual void addEdge(OldElement* edge)
	{
		_edges.push_back(edge);
	}

protected:
	void setFace(size_t index, OldElement* face) { _faces[index] = face; }
	void setEdge(size_t index, OldElement* edge) { _edges[index] = edge; }

	eslocal _params[PARAMS_SIZE];
	std::vector<OldElement*> _edges;
	std::vector<OldElement*> _faces;
};

}

#endif /* SRC_MESH_ELEMENTS_VOLUME_VOLUMEELEMENT_H_ */
