
#include "monitoring.h"

#include "../../../config/ecf/output.h"
#include "../../../assembler/step.h"
#include "../../../mesh/mesh.h"

//#include "../../config/ecf/output.h"
//#include "../../old/mesh/settings/property.h"
//#include "../../old/mesh/structures/mesh.h"
//#include "../../old/mesh/structures/region.h"
//#include "../../assembler/step.h"
//#include "../../assembler/solution.h"
//
//#include "../../config/ecf/environment.h"
//#include "../../basis/logging/logging.h"
//#include "../../basis/utilities/parser.h"
//
//#include <cmath>

using namespace espreso;

char Monitoring::delimiter = ';';

std::string center(const std::string &value, size_t size)
{
	std::string ret;
	ret.insert(ret.end(), (size - value.size()) / 2 + (size - value.size()) % 2, ' ');
	ret.insert(ret.end(), value.begin(), value.end());
	ret.insert(ret.end(), (size - value.size()) / 2, ' ');
	return ret;
}

std::string right(const std::string &value, size_t size)
{
	std::string ret;
	ret.insert(ret.end(), size - value.size(), ' ');
	ret.insert(ret.end(), value.begin(), value.end());
	return ret;
}


Monitoring::Monitoring(const Mesh &mesh, const OutputConfiguration &configuration)
: ResultStoreBase(mesh), _mesh(mesh)
{

}

void Monitoring::updateSolution(const Step &step)
{

}

//std::vector<espreso::Property> Monitoring::getProperties(const std::string &name)
//{
//	for (auto p = _mesh->propertyGroups().begin(); p != _mesh->propertyGroups().end(); ++p) {
//		if (p->second.size() > 1) {
//			std::stringstream ss; ss << p->first;
//			std::string pname = ss.str().substr(0, ss.str().find_last_of("_"));
//			if (StringCompare::caseInsensitiveEq(pname, name)) {
//				return p->second;
//			}
//		}
//	}
//
//	std::stringstream ss(name);
//	espreso::Property property;
//	ss >> property;
//
//	return { property };
//}
//
//Monitor::Monitor()
//: printSize(10), region(NULL), statistics(StatisticalData::EMPTY)
//{
//
//}

//void Monitoring::updateMesh()
//{
//	_monitors.reserve(_configuration.monitoring.size());
//	for (auto it = _configuration.monitoring.begin(); it != _configuration.monitoring.end(); ++it) {
//		_monitors.resize(it->first);
//		_monitors.back().region = _mesh->region(it->second.region);
//		_mesh->addMonitoredRegion(_monitors.back().region);
//		_monitors.back().statistics = static_cast<StatisticalData>(std::pow(2, static_cast<int>(it->second.statistics)));
//		_monitors.back().properties = getProperties(it->second.property);
//		_monitors.back().printSize = std::max(std::max((size_t)10, it->second.region.size()), it->second.property.size()) + 4;
//	}
//
//	if (environment->MPIrank) {
//		return;
//	}
//
//	_os.open(Logging::outputRoot() + "/" + Logging::name + ".emr");
//	if (!_os.is_open()) {
//		ESINFO(GLOBAL_ERROR) << "Cannot open file for storing monitor report\n";
//	}
//
//	_os << "\n";
//	_os << std::string(9, ' ') << delimiter << std::string(9, ' ') << delimiter;
//	for (size_t i = 0; i < _monitors.size(); i++) {
//		if (_monitors[i].region == NULL) {
//			_os << center("---", _monitors[i].printSize) << delimiter;
//		} else {
//			_os << center(_monitors[i].region->name, _monitors[i].printSize) << delimiter;
//		}
//	}
//	_os << "\n";
//
//	_os << right("step", 9) << delimiter << right("substep", 9) << delimiter;
//	for (size_t i = 0; i < _monitors.size(); i++) {
//		std::stringstream ss;
//		if (_monitors[i].properties.size() == 0) {
//			ss << "-";
//		} else if (_monitors[i].properties.size() > 1) {
//			std::stringstream ssp; ssp << _monitors[i].properties[0];
//			ss << ssp.str().substr(0, ssp.str().find_last_of("_"));
//		} else {
//			ss << _monitors[i].properties[0];
//		}
//		_os << center(ss.str(), _monitors[i].printSize) << delimiter;
//	}
//	_os << "\n";
//
//	_os << std::string(9, ' ') << delimiter << std::string(9, ' ') << delimiter;
//	for (size_t i = 0; i < _monitors.size(); i++) {
//		switch (_monitors[i].statistics) {
//		case StatisticalData::EMPTY:   _os << center("<--->"    , _monitors[i].printSize) << delimiter; break;
//		case StatisticalData::AVERAGE: _os << center("<AVERAGE>", _monitors[i].printSize) << delimiter; break;
//		case StatisticalData::MIN:     _os << center("<MIN>"    , _monitors[i].printSize) << delimiter; break;
//		case StatisticalData::MAX:     _os << center("<MAX>"    , _monitors[i].printSize) << delimiter; break;
//		case StatisticalData::NORM:    _os << center("<NORM>"   , _monitors[i].printSize) << delimiter; break;
//		default: break;
//		}
//	}
//	_os << "\n\n";
//
//}
//void Monitoring::storeSolution(const Step &step, const std::vector<Solution*> &solution, const std::vector<std::pair<ElementType, Property> > &properties)
//{
//	switch (_configuration.monitors_store_frequency) {
//	case OutputConfiguration::STORE_FREQUENCY::NEVER:
//		return;
//
//	case OutputConfiguration::STORE_FREQUENCY::EVERY_NTH_TIMESTEP:
//		if ((step.substep + 1) % _configuration.monitors_nth_stepping != 0) {
//			return;
//		}
//		break;
//	case OutputConfiguration::STORE_FREQUENCY::LAST_TIMESTEP:
//		if (!step.isLast()) {
//			return;
//		}
//		break;
//	case OutputConfiguration::STORE_FREQUENCY::EVERY_TIMESTEP:
//	case OutputConfiguration::STORE_FREQUENCY::DEBUG:
//		break;
//	}
//
//	for (size_t i = 0; i < _monitors.size(); i++) {
//		if (_monitors[i].statistics == StatisticalData::EMPTY) {
//			continue;
//		}
//		for (size_t s = 0; s < solution.size(); s++) {
//			if (solution[s]->hasProperty(_monitors[i].properties[0])) {
//				solution[s]->computeStatisticalData(step);
//			}
//		}
//	}
//
//	if (environment->MPIrank) {
//		return;
//	}
//	_os << right(std::to_string(step.step + 1), 8) << " " << delimiter;
//	_os << right(std::to_string(step.substep + 1), 8) << " " << delimiter;
//
//	for (size_t i = 0; i < _monitors.size(); i++) {
//		if (_monitors[i].statistics == StatisticalData::EMPTY) {
//			_os << center("|", _monitors[i].printSize) << delimiter;
//			continue;
//		}
//		double value;
//		bool found = false;
//		for (size_t s = 0; s < solution.size(); s++) {
//			if (solution[s]->hasProperty(_monitors[i].properties[0])) {
//				value = solution[s]->getStatisticalData(_monitors[i].properties, _monitors[i].statistics, _monitors[i].region);
//				found = true;
//				break;
//			}
//		}
//		if (!found) {
//			ESINFO(GLOBAL_ERROR) << "ESPRESO monitor error: request for unknown property: " << _monitors[i].properties[0];
//		}
//
//		std::stringstream ss;
//		ss << std::scientific << value;
//		_os << center(ss.str(), _monitors[i].printSize) << delimiter;
//	}
//	_os << "\n";
//
//}
//
//void Monitoring::finalize()
//{
//	if (environment->MPIrank) {
//		return;
//	}
//	_os.close();
//}







