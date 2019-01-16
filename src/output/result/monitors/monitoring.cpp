

#include "monitoring.h"

#include "esinfo/time.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/mpiinfo.h"

#include "basis/containers/serializededata.h"
#include "basis/logging/logging.h"
#include "basis/utilities/utils.h"
#include "basis/utilities/parser.h"

#include "config/ecf/output.h"

#include "mesh/mesh.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/boundaryregionstore.h"
#include "mesh/store/elementsregionstore.h"

#include <iomanip>

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

bool Monitoring::storeStep()
{
	switch (info::ecf->output.monitors_store_frequency) {
	case OutputConfiguration::STORE_FREQUENCY::NEVER:
		return false;
	case OutputConfiguration::STORE_FREQUENCY::EVERY_TIMESTEP:
		return true;
	case OutputConfiguration::STORE_FREQUENCY::EVERY_NTH_TIMESTEP:
		return time::substep % info::ecf->output.monitors_nth_stepping == 0;
	case OutputConfiguration::STORE_FREQUENCY::LAST_TIMESTEP:
		return time::isLast();
	default:
		return false;
	}
}


Monitoring::Monitoring(const Mesh &mesh)
: ResultStoreBase(mesh)
{

}

Monitoring::~Monitoring()
{

}

void Monitoring::updateMesh()
{
	for (auto it = info::ecf->output.monitoring.begin(); it != info::ecf->output.monitoring.end(); ++it) {
		if (it->first <= 0) {
			ESINFO(GLOBAL_ERROR) << "Invalid column index in monitoring.";
		}
		if (it->first > _monitors.size()) {
			_monitors.resize(it->first);
		}

		ElementsRegionStore *estore = NULL;
		BoundaryRegionStore *bstore = NULL;
		bool regionNotFound = true;
		for (size_t r = 0; regionNotFound && r < _mesh.elementsRegions.size(); r++) {
			if (StringCompare::caseInsensitiveEq(it->second.region, _mesh.elementsRegions[r]->name)) {
				estore = _mesh.elementsRegions[r];
				regionNotFound = false;
			}
		}
		for (size_t r = 0; regionNotFound && r < _mesh.boundaryRegions.size(); r++) {
			if (StringCompare::caseInsensitiveEq(it->second.region, _mesh.boundaryRegions[r]->name)) {
				bstore = _mesh.boundaryRegions[r];
				regionNotFound = false;
			}
		}
		if (regionNotFound) {
			ESINFO(GLOBAL_ERROR) << "Monitoring contains unknown region '" << it->second.region << "'.";
		}

		NodeData *ndata = NULL;
		ElementData *edata = NULL;
		bool propertyNotFound = true;
		for (size_t i = 0; propertyNotFound && i < _mesh.nodes->data.size(); i++) {
			for (size_t p = 0; p < _mesh.nodes->data[i]->names.size(); p++) {
				if (StringCompare::caseInsensitiveEq(it->second.property, _mesh.nodes->data[i]->names[p])) {
					ndata = _mesh.nodes->data[i];
					propertyNotFound = false;
					break;
				}
			}
		}

		for (size_t i = 0; propertyNotFound && i < _mesh.elements->data.size(); i++) {
			for (size_t p = 0; p < _mesh.elements->data[i]->names.size(); p++) {
				if (StringCompare::caseInsensitiveEq(it->second.property, _mesh.elements->data[i]->names[p])) {
					edata = _mesh.elements->data[i];
					propertyNotFound = false;
					break;
				}
			}
		}
		if (propertyNotFound) {
			ESINFO(GLOBAL_ERROR) << "Monitoring contains unknown property '" << it->second.property << "'.";
		}

		if (edata != NULL && bstore != NULL) {
			ESINFO(GLOBAL_ERROR) << "Cannot monitor element property '" << it->second.property << "' on element region '" << it->second.region << "'.";
		}
		if (edata != NULL && estore != NULL) {
			for (size_t i = 0; i < _edata.size(); i++) {
				if (_edata[i].first->names == edata->names && _edata[i].second->name == estore->name) {
					edata = NULL;
					estore = NULL;
					break;
				}
			}
		}
		if (ndata != NULL && bstore != NULL) {
			for (size_t i = 0; i < _nbdata.size(); i++) {
				if (_nbdata[i].first->names == ndata->names && _nbdata[i].second->name == bstore->name) {
					ndata = NULL;
					bstore = NULL;
					break;
				}
			}
		}
		if (ndata != NULL && estore != NULL) {
			for (size_t i = 0; i < _nedata.size(); i++) {
				if (_nedata[i].first->names == ndata->names && _nedata[i].second->name == estore->name) {
					ndata = NULL;
					estore = NULL;
					break;
				}
			}
		}

		if (edata != NULL && estore != NULL) {
			_edata.push_back(std::make_pair(edata, estore));
		}
		if (ndata != NULL && bstore != NULL) {
			_nbdata.push_back(std::make_pair(ndata, bstore));
		}
		if (ndata != NULL && estore != NULL) {
			_nedata.push_back(std::make_pair(ndata, estore));
		}
	}

	for (size_t i = 0; i < _edata.size(); i++) {
		_statistics.resize(_statistics.size() + _edata[i].first->names.size());
	}
	for (size_t i = 0; i < _nbdata.size(); i++) {
		_statistics.resize(_statistics.size() + _nbdata[i].first->names.size());
	}
	for (size_t i = 0; i < _nedata.size(); i++) {
		_statistics.resize(_statistics.size() + _nedata[i].first->names.size());
	}

	for (auto it = info::ecf->output.monitoring.begin(); it != info::ecf->output.monitoring.end(); ++it) {
		_monitors[it->first - 1].name = it->second.region;
		_monitors[it->first - 1].property = it->second.property;
		_monitors[it->first - 1].printSize = std::max(std::max((size_t)10, it->second.region.size()), it->second.property.size()) + 4;

		esint offset = 0;
		for (size_t i = 0; i < _edata.size(); i++) {
			for (size_t p = 0; p < _edata[i].first->names.size(); ++p, ++offset) {
				if (
						StringCompare::caseInsensitiveEq(it->second.property, _edata[i].first->names[p]) &&
						StringCompare::caseInsensitiveEq(it->second.region, _edata[i].second->name)) {

					_monitors[it->first - 1].data = (double*)(_statistics.data() + offset);
				}
			}
		}
		for (size_t i = 0; i < _nbdata.size(); i++) {
			for (size_t p = 0; p < _nbdata[i].first->names.size(); ++p, ++offset) {
				if (
						StringCompare::caseInsensitiveEq(it->second.property, _nbdata[i].first->names[p]) &&
						StringCompare::caseInsensitiveEq(it->second.region, _nbdata[i].second->name)) {

					_monitors[it->first - 1].data = (double*)(_statistics.data() + offset);
				}
			}
		}
		for (size_t i = 0; i < _nedata.size(); i++) {
			for (size_t p = 0; p < _nedata[i].first->names.size(); ++p, ++offset) {
				if (
						StringCompare::caseInsensitiveEq(it->second.property, _nedata[i].first->names[p]) &&
						StringCompare::caseInsensitiveEq(it->second.region, _nedata[i].second->name)) {

					_monitors[it->first - 1].data = (double*)(_statistics.data() + offset);
				}
			}
		}

		switch (it->second.statistics) {
		case MonitorConfiguration::STATISTICS::MIN:
			_monitors[it->first - 1].stats = "<MIN>"; break;
		case MonitorConfiguration::STATISTICS::MAX:
			_monitors[it->first - 1].data += 1;
			_monitors[it->first - 1].stats = "<MAX>"; break;
		case MonitorConfiguration::STATISTICS::AVG:
			_monitors[it->first - 1].data += 2;
			_monitors[it->first - 1].stats = "<AVERAGE>"; break;
		case MonitorConfiguration::STATISTICS::NORM:
			_monitors[it->first - 1].data += 3;
			_monitors[it->first - 1].stats = "<NORM>"; break;
		}
	}

	_os.open(Logging::outputRoot() + "/" + Logging::name + ".emr");

	if (!_os.is_open()) {
		ESINFO(GLOBAL_ERROR) << "Cannot open file for storing monitor report\n";
	}

	_os << "\n";
	_os << std::string(9, ' ') << delimiter << std::string(9, ' ') << delimiter << std::string(9, ' ') << delimiter;
	for (size_t i = 0; i < _monitors.size(); i++) {
		_os << center(_monitors[i].name, _monitors[i].printSize) << delimiter;
	}
	_os << "\n";

	_os << right("step ", 9) << delimiter << right("substep ", 9) << delimiter << right("time ", 9) << delimiter;
	for (size_t i = 0; i < _monitors.size(); i++) {
		_os << center(_monitors[i].property, _monitors[i].printSize) << delimiter;
	}
	_os << "\n";

	_os << std::string(9, ' ') << delimiter << std::string(9, ' ') << delimiter << std::string(9, ' ') << delimiter;
	for (size_t i = 0; i < _monitors.size(); i++) {
		_os << center(_monitors[i].stats, _monitors[i].printSize) << delimiter;
	}
	_os << "\n\n";
	_os.flush();
}

void Monitoring::updateSolution()
{
	if (!storeStep()) {
		return;
	}

	esint offset = 0;
	for (size_t i = 0; i < _edata.size(); offset += _edata[i++].first->names.size()) {
		_edata[i].first->statistics(_edata[i].second->elements->datatarray(), _edata[i].second->uniqueTotalSize, _statistics.data() + offset);
	}
	for (size_t i = 0; i < _nbdata.size(); offset += _nbdata[i++].first->names.size()) {
		_nbdata[i].first->statistics(_nbdata[i].second->nodes->datatarray(), _nbdata[i].second->uniqueTotalSize, _statistics.data() + offset);
	}
	for (size_t i = 0; i < _nedata.size(); offset += _nedata[i++].first->names.size()) {
		_nedata[i].first->statistics(_nedata[i].second->nodes->datatarray(), _nedata[i].second->uniqueTotalSize, _statistics.data() + offset);
	}

	if (info::mpi::MPIrank) {
		return;
	}

	_os << right(std::to_string(time::step + 1), 8) << " " << delimiter;
	_os << right(std::to_string(time::substep + 1), 8) << " " << delimiter;
	std::stringstream time; time << std::setprecision(4) << std::fixed << time::current;
	_os << right(time.str(), 8) << " " << delimiter;

	for (size_t i = 0; i < _monitors.size(); i++) {
		std::stringstream value;
		if (_monitors[i].data != NULL) {
			value << std::scientific << *_monitors[i].data;
		}
		_os << center(value.str(), _monitors[i].printSize) << delimiter;
	}
	_os << "\n";
	_os.flush();
}







