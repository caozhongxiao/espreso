
#include "ensight.h"

#include "../../../basis/containers/point.h"
#include "../../../basis/containers/serializededata.h"
#include "../../../basis/logging/logging.h"

#include "../../../mesh/mesh.h"
#include "../../../mesh/store/nodestore.h"
#include "../../../mesh/store/elementstore.h"

#include <iomanip>

using namespace espreso;

EnSight::EnSight(const std::string &name, const Mesh &mesh)
: _path(Logging::outputRoot() + "/"), _name(name), _mesh(mesh),
  _casefile(_path + _name + ".case")
{
	_casefile << "FORMAT\n";
	_casefile << "type: ensight gold\n";

	_casefile << "GEOMETRY\n";
	_casefile << "model:\t" << _name << ".geo\n";

//	_os << "VARIABLE\n";
	_casefile.flush();
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
	os << "COORDINATES (x1, x2, x3, ..., y1, y2, y3, ..., z1, z2, z3, ...)\n";

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
//	for (size_t r = 0; r < _meshInfo->regions(); r++) {
//		for (size_t n = 0; n < _meshInfo->region(r).data.pointDataDouble.size(); n++) {
//
//		}
//	}
//
//	std::ofstream os(name + ".temperature");
//	os << "TEMPERATURE\n";
//	for (size_t r = 0; r < _meshInfo->regions(); r++) {
//		os << "part\n";
//		os << std::setw(10) << r + 1 << "\n";
//
//		os << "coordinates\n";
//		for (size_t n = 0; n < _meshInfo->region(r).data.pointDataDouble.size(); n++) {
//			os << std::scientific << std::setprecision(5) << " " << _meshInfo->region(r).coordinates[n] << "\n";
//		}
//
//		os << "hexa8\n";
//		os << std::setw(10) << _meshInfo->region(r).elementsTypes.size() << "\n";
//		for (size_t e = 0; e < _meshInfo->region(r).elementsTypes.size(); e++) {
//			for (size_t n = 0; n < 8; n++) {
//				os << std::setw(10) << _meshInfo->region(r).elements[8 * e + n] + 1;
//			}
//			os << "\n";
//		}
//	}
//	os.close();
}




