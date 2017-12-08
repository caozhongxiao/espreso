
#include "solution.h"

#include "../basis/logging/logging.h"

#include "../mesh/mesh.h"
#include "../mesh/store/elementstore.h"
#include "../mesh/store/nodestore.h"

// TODO: MESH
#include "../old/mesh/structures/elementtypes.h"

using namespace espreso;

Solution::Solution(const Mesh &mesh, const std::string &name, ElementType eType, size_t DOFs, std::vector<std::vector<double> > &data)
: name(name), eType(eType), DOFs(DOFs), data(data)
{

}

Solution::Solution(const Mesh &mesh, const std::string &name, ElementType eType, size_t DOFs)
: name(name), eType(eType), DOFs(DOFs), data(_data)
{
	// TODO: MESH
	_data.resize(mesh.elements->ndomains);
	for (size_t p = 0; p < mesh.elements->ndomains; p++) {
		switch (eType) {
		case ElementType::ELEMENTS:
			_data[p].resize(DOFs * (mesh.elements->elementsDistribution[p + 1] - mesh.elements->elementsDistribution[p]));
			break;
		case ElementType::NODES:
			_data[p].resize(DOFs * (mesh.nodes->dintervals[p].back().DOFOffset + mesh.nodes->dintervals[p].back().end - mesh.nodes->dintervals[p].back().begin));
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: invalid solution element type.";
		}
	}
}

void Solution::fill(double value)
{
	for (size_t p = 0; p < _data.size(); p++) {
		std::fill(_data[p].begin(), _data[p].end(), value);
	}
}


