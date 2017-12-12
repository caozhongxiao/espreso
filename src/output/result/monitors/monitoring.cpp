
#include "monitoring.h"

#include "../../../basis/logging/logging.h"
#include "../../../basis/utilities/utils.h"
#include "../../../basis/utilities/parser.h"

#include "../../../config/ecf/environment.h"
#include "../../../config/ecf/output.h"
#include "../../../assembler/step.h"

#include "../../../mesh/mesh.h"
#include "../../../mesh/store/nodestore.h"
#include "../../../mesh/store/boundaryregionstore.h"
#include "../../../mesh/store/elementsregionstore.h"

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


Monitoring::Monitoring(const Mesh &mesh, const OutputConfiguration &configuration, bool async)
: ResultStoreBase(mesh), _configuration(configuration), _async(async)
{
	MPI_Comm_split(environment->MPICommunicator, 0, environment->MPIrank, &_communicator);
}

Monitoring::~Monitoring()
{
	MPI_Comm_free(&_communicator);
}

void Monitoring::updateMesh()
{
	for (auto it = _configuration.monitoring.begin(); it != _configuration.monitoring.end(); ++it) {
		if (it->first <= 0) {
			ESINFO(GLOBAL_ERROR) << "Invalid column index in monitoring.";
		}
		if (it->first > _monitors.size()) {
			_monitors.resize(it->first);
		}
		for (size_t r = 0; r < _mesh.elementsRegions.size(); r++) {
			if (StringCompare::caseInsensitiveEq(it->second.region, _mesh.elementsRegions[r]->name)) {
				_eregions.push_back(_mesh.elementsRegions[r]);
			}
		}
		for (size_t r = 0; r < _mesh.boundaryRegions.size(); r++) {
			if (StringCompare::caseInsensitiveEq(it->second.region, _mesh.boundaryRegions[r]->name)) {
				_bregions.push_back(_mesh.boundaryRegions[r]);
			}
		}
	}
	Esutils::sortAndRemoveDuplicity(_eregions);
	Esutils::sortAndRemoveDuplicity(_bregions);

	_data.resize(_eregions.size() + _bregions.size());

	for (auto it = _configuration.monitoring.begin(); it != _configuration.monitoring.end(); ++it) {
		_monitors[it->first - 1].name = it->second.region;
		_monitors[it->first - 1].property = it->second.property;
		_monitors[it->first - 1].printSize = std::max(std::max((size_t)10, it->second.region.size()), it->second.property.size()) + 4;
		for (size_t r = 0; r < _eregions.size(); r++) {
			if (StringCompare::caseInsensitiveEq(it->second.region, _eregions[r]->name)) {
				switch (it->second.statistics) {
				case MonitorConfiguration::STATISTICS::MIN:
					_monitors[it->first - 1].stats = "<MIN>";
					_monitors[it->first - 1].data = &_data[r].min; continue;
				case MonitorConfiguration::STATISTICS::MAX:
					_monitors[it->first - 1].stats = "<MAX>";
					_monitors[it->first - 1].data = &_data[r].max; continue;
				case MonitorConfiguration::STATISTICS::AVG:
					_monitors[it->first - 1].stats = "<AVERAGE>";
					_monitors[it->first - 1].data = &_data[r].avg; continue;
				case MonitorConfiguration::STATISTICS::NORM:
					_monitors[it->first - 1].stats = "<NORM>";
					_monitors[it->first - 1].data = &_data[r].norm; continue;
				}
			}
		}
		for (size_t r = 0; r < _bregions.size(); r++) {
			if (StringCompare::caseInsensitiveEq(it->second.region, _bregions[r]->name)) {
				switch (it->second.statistics) {
				case MonitorConfiguration::STATISTICS::MIN:
					_monitors[it->first - 1].stats = "<MIN>";
					_monitors[it->first - 1].data = &_data[_eregions.size() + r].min; continue;
				case MonitorConfiguration::STATISTICS::MAX:
					_monitors[it->first - 1].stats = "<MAX>";
					_monitors[it->first - 1].data = &_data[_eregions.size() + r].max; continue;
				case MonitorConfiguration::STATISTICS::AVG:
					_monitors[it->first - 1].stats = "<AVERAGE>";
					_monitors[it->first - 1].data = &_data[_eregions.size() + r].avg; continue;
				case MonitorConfiguration::STATISTICS::NORM:
					_monitors[it->first - 1].stats = "<NORM>";
					_monitors[it->first - 1].data = &_data[_eregions.size() + r].norm; continue;
				}
			}
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

void Monitoring::updateSolution(const Step &step)
{
	eslocal offset = 0;
	for (size_t i = 0; i < _eregions.size(); ++i, ++offset) {
		_mesh.computeGatheredNodeStatistic(_mesh.nodes->data.front(), _eregions[i], _data.data() + offset, _communicator);
	}
	for (size_t i = 0; i < _bregions.size(); ++i, ++offset) {
		_mesh.computeGatheredNodeStatistic(_mesh.nodes->data.front(), _bregions[i], _data.data() + offset, _communicator);
	}

	if (environment->MPIrank) {
		return;
	}

	_os << right(std::to_string(step.step + 1), 8) << " " << delimiter;
	_os << right(std::to_string(step.substep + 1), 8) << " " << delimiter;
	std::stringstream time; time << std::setprecision(4) << std::fixed << step.currentTime;
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







