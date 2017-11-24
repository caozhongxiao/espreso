
#include "ensight.h"

#include "../../../basis/containers/point.h"
#include "../../../basis/containers/serializededata.h"
#include "../../../basis/logging/logging.h"

#include "../../../config/ecf/environment.h"

#include "../../../mesh/mesh.h"
#include "../../../mesh/store/nodestore.h"
#include "../../../mesh/store/elementstore.h"

#include <fstream>
#include <iomanip>

using namespace espreso;

EnSight::EnSight(const std::string &name, const Mesh &mesh)
: _path(Logging::outputRoot() + "/"), _name(name), _mesh(mesh), _casefile(NULL)
{
	if (environment->MPIrank == 0) {
		_casefile = new std::ofstream(_path + _name + ".case");
		(*_casefile) << "#\n";
		(*_casefile) << "# ESPRESO solution\n";
		(*_casefile) << "#\n";

		(*_casefile) << "\nFORMAT\n";
		(*_casefile) << "type: \tensight gold\n\n";

		(*_casefile) << "GEOMETRY\n\n";
		(*_casefile) << "model:\t" << _name << ".geo\n\n";

		(*_casefile).flush();
	}
}

EnSight::~EnSight()
{
	if (_casefile != NULL) {
		delete _casefile;
	}
}

void EnSight::storeGeometry()
{
	std::ofstream os(_path + _name + ".geo");
	os << "EnSight Gold geometry format\n";
	os << "----------------------------\n";

	os << "node id off\n";
	os << "element id off\n";

	os << "part\n";
	os << std::setw(10) << 1 << "\n";
	os << "MESH\n";

	os << "coordinates\n";
	os << std::setw(10) << _mesh.nodes->coordinates->datatarray().size() << "\n";

	for (auto n = _mesh.nodes->coordinates->datatarray().begin(); n != _mesh.nodes->coordinates->datatarray().end(); ++n) {
		os << std::scientific << std::setprecision(5) << " " << n->x << "\n";
	}
	for (auto n = _mesh.nodes->coordinates->datatarray().begin(); n != _mesh.nodes->coordinates->datatarray().end(); ++n) {
		os << std::scientific << std::setprecision(5) << " " << n->y << "\n";
	}
	for (auto n = _mesh.nodes->coordinates->datatarray().begin(); n != _mesh.nodes->coordinates->datatarray().end(); ++n) {
		os << std::scientific << std::setprecision(5) << " " << n->z << "\n";
	}

	os << "hexa8\n";
	os << std::setw(10) << _mesh.elements->size << "\n";
	auto node = _mesh.elements->nodes->datatarray().cbegin();
	for (size_t e = 0; e < _mesh.elements->size; e++) {
		for (size_t n = 0; n < 8; ++n, ++node) {
			os << std::setw(10) << *node + 1;
		}
		os << "\n";
	}
	os.close();
}

void EnSight::storeVariables()
{
	if (environment->MPIrank == 0) {
		(*_casefile) << "VARIABLE\n\n";
		for (size_t i = 0; i < _mesh.nodes->data.size(); i++) {
			std::string filename = _name + "." + _mesh.nodes->data[i]->names.front().substr(0, 4);
			(*_casefile) << "scalar per node:\t" << _mesh.nodes->data[i]->names.front() << "\t" << filename << "\n";
		}
		_casefile->flush();
	}
	for (size_t i = 0; i < _mesh.nodes->data.size(); i++) {
		std::string filename = _name + "." + _mesh.nodes->data[i]->names.front().substr(0, 4);

		std::ofstream os(_path + filename);
		os << _mesh.nodes->data[i]->names.front() << "\n";
		os << "part\n";
		os << std::setw(10) << 1 << "\n";
		os << "coordinates\n";

		for (size_t n = 0; n < _mesh.nodes->data[i]->gathredData->size(); ++n) {
			os << std::scientific << std::setprecision(5) << " " << (*_mesh.nodes->data[i]->gathredData)[n] << "\n";
		}

		os.close();
	}
}




