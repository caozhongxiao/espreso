
#include "output.h"

#include "newmesh.h"
#include "elements/newelement.h"
#include "elements/elementstore.h"

#include "../basis/point/point.h"
#include "../basis/containers/serializededata.h"

#include "../config/ecf/environment.h"

#include <fstream>

using namespace espreso;

void NewOutput::VTKLegacy(const std::string &name, ElementStore *elements, ElementStore *nodes)
{
	std::ofstream os(name + std::to_string(environment->MPIrank) + ".vtk");

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
	os << "\n";

	os << "CELL_TYPES " << elements->size << "\n";
	for (auto e = elements->epointers->datatarray().begin(); e != elements->epointers->datatarray().end(); ++e) {
		switch ((*e)->code) {
		case NewElement::CODE::SQUARE4:
			os << "9\n";
			break;
		case NewElement::CODE::SQUARE8:
			os << "23\n";
			break;
		case NewElement::CODE::TRIANGLE3:
			os << "5\n";
			break;
		case NewElement::CODE::TRIANGLE6:
			os << "22\n";
			break;
		case NewElement::CODE::TETRA4:
			os << "10\n";
			break;
		case NewElement::CODE::TETRA10:
			os << "24\n";
			break;
		case NewElement::CODE::PYRAMID5:
			os << "14\n";
			break;
		case NewElement::CODE::PYRAMID13:
			os << "27\n";
			break;
		case NewElement::CODE::PRISMA6:
			os << "13\n";
			break;
		case NewElement::CODE::PRISMA15:
			os << "26\n";
			break;
		case NewElement::CODE::HEXA8:
			os << "12\n";
			break;
		case NewElement::CODE::HEXA20:
			os << "25\n";
			break;
		}
	}
	os << "\n";

	os << "CELL_DATA " << elements->size << "\n";
	os << "SCALARS cluster int 1\n";
	os << "LOOKUP_TABLE default\n";
	for (size_t e = 0; e < elements->size; e++) {
		os << environment->MPIrank << "\n";
	}
}



