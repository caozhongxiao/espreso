
#include "output.h"

#include "newmesh.h"
#include "elements/elementstore.h"

#include "../basis/point/point.h"
#include "../basis/containers/serializededata.h"

#include "../config/ecf/environment.h"

#include <fstream>

using namespace espreso;

void NewOutput::VTKLegacy(ElementStore *elements, ElementStore *nodes)
{
	std::ofstream os("res" + std::to_string(environment->MPIrank) + ".vtk");

	os << "# vtk DataFile Version 2.0\n";
	os << "EXAMPLE\n";
	os << "ASCII\n";
	os << "DATASET UNSTRUCTURED_GRID\n\n";

	os << "POINTS " << nodes->size << " float\n";
	for (auto n = nodes->coordinates->datatarray().begin(); n != nodes->coordinates->datatarray().end(); ++n) {
		os << n->x << " " << n->y << " " << n->z << "\n";
	}
	os << "\n";

	os << "CELLS " << elements->size << " " << elements->size + elements->nodes->datatarray().size() << "\n";
	for (auto e = elements->nodes->cbegin(); e != elements->nodes->cend(); ++e) {
		os << e->size() << " ";
		for (auto n = e->begin(); n != e->end(); ++n) {
			os << *n << " ";
		}
		os << "\n";
	}
}



