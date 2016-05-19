
#include "configuration.h"

using namespace espreso::input;

static void fillParameters(std::map<std::string, Parameter*, CaseInsensitiveCompare> &parameters, const std::vector<Description> &description)
{
	for (size_t i = 0; i < description.size(); i++) {
		switch (description[i].type) {
		case STRING_PARAMETER: {
			parameters[description[i].name] = new StringParameter(
				description[i].name,
				description[i].description
			);
			break;
		}
		case INTEGER_PARAMETER: {
			parameters[description[i].name] = new IntegerParameter(
				description[i].name,
				description[i].description
			);
			break;
		}
		case DOUBLE_PARAMETER: {
			parameters[description[i].name] = new DoubleParameter(
				description[i].name,
				description[i].description
			);
			break;
		}
		case BOOLEAN_PARAMETER: {
			parameters[description[i].name] = new BooleanParameter(
				description[i].name,
				description[i].description
			);
			break;
		}
		}
	}
}

Configuration::Configuration(std::vector<Description> &description, const Options &options)
{
	fillParameters(_parameters, description);

	_parameters["CMD_LINE_ARGUMENTS"] = new StringParameter(
		"CMD_LINE_ARGUMENTS",
		"Arguments set from the command line."
	);

	load(options);
	set(description);
}

Configuration::Configuration(std::vector<Description> &description, const std::string &path)
{
	std::ifstream file(path);

	if (file.is_open()) {
		// check existence of the file
		file.close();

		fillParameters(_parameters, description);
		Options options;
		options.path = path;

		load(options);
		set(description);
	}
}

void Configuration::load(const Options &options)
{
	std::ifstream file(options.path);
	std::map<std::string, Parameter*>::iterator it;

	if (file.is_open()) {
		std::string line;

		while (getline(file, line, '\n')) {
			line.erase(0, line.find_first_not_of(" "));
			if (!line.size() || line.compare(0, 1, "#") == 0) {
				continue;
			}
			size_t pos = line.find("#");
			if (pos < std::string::npos) {
				line = line.substr(0, pos);
			}
			for (it = _parameters.begin(); it != _parameters.end(); ++it) {;
				if (it->second->match(line)) {
					if (it->second->isSet()) {
						ESINFO(ALWAYS) << "Warning: parameter " << it->second->name() << " is set more than once.";
					}
					it->second->set(line);
					break;
				}
			}
		}

		file.close();
	} else {
		ESINFO(ERROR) << "The example on path '" << options.path << "' not found.";
	}

	// Read attributes from command line
	if (_parameters.find("CMD_LINE_ARGUMENTS") == _parameters.end()) {
		return;
	}

	std::vector<std::pair<std::string, int> > cmdLine;
	size_t cmdLineSize = 0;
	if (_parameters["CMD_LINE_ARGUMENTS"]->isSet()) {
		std::string val = value<std::string>("CMD_LINE_ARGUMENTS", "");
		while(true) {
			size_t pos = val.find(" ");
			std::string argument = val.substr(0, pos);
			for (it = _parameters.begin(); it != _parameters.end(); ++it) {
				if (it->second->match(argument)) {
					cmdLine.push_back(std::pair<std::string, int>(it->second->name(), cmdLineSize));
					break;
				}
			}
			cmdLineSize++;
			val = val.erase(0, pos);
			val = val.erase(0, val.find_first_not_of(" "));
			if (!val.size()) {
				break;
			}
		}
	}

	if (options.nameless.size() < cmdLineSize) {
		ESINFO(ERROR) << "Too few command line arguments. ESPRESO assumes " << value<std::string>("CMD_LINE_ARGUMENTS", "") << "\n";
		exit(EXIT_FAILURE);
	}
	if (options.nameless.size() > cmdLineSize) {
		ESINFO(ALWAYS) << "Warning: ESPRESO omits some command line arguments.\n";
	}
	for (size_t i = 0; i < cmdLine.size(); i++) {
		_parameters[cmdLine[i].first]->set(std::string(_parameters[cmdLine[i].first]->name() + "=" + options.nameless[cmdLine[i].second]));
	}
}

void Configuration::set(std::vector<Description> &description)
{
	for (size_t i = 0; i < description.size(); i++) {
		auto it = _parameters.find(description[i].name);
		if (it != _parameters.end() && it->second->isSet()) {
			switch (description[i].type) {
			case INTEGER_PARAMETER:
				*(static_cast<int*>(description[i].value)) = static_cast<IntegerParameter*>(it->second)->get();
				break;
			case LONG_PARAMETER:
				*(static_cast<long*>(description[i].value)) = static_cast<IntegerParameter*>(it->second)->get();
				break;
			case SIZE_PARAMETER:
				*(static_cast<size_t*>(description[i].value)) = static_cast<IntegerParameter*>(it->second)->get();
				break;
			case DOUBLE_PARAMETER:
				*(static_cast<double*>(description[i].value)) = static_cast<DoubleParameter*>(it->second)->get();
				break;
			case STRING_PARAMETER:
				*(static_cast<std::string*>(description[i].value)) = static_cast<StringParameter*>(it->second)->get();
				break;
			case BOOLEAN_PARAMETER:
				*(static_cast<bool*>(description[i].value)) = static_cast<BooleanParameter*>(it->second)->get();
				break;
			}
		}
	}

}

void Configuration::print() const
{
	std::map<std::string, Parameter*>::const_iterator it;
	for (it = _parameters.begin(); it != _parameters.end(); ++it) {
		if (!it->second->isSet()) {
			continue;
		}
		std::stringstream ss;
		ss << it->second->name() << " = ";
		switch (it->second->type()) {
			case STRING_PARAMETER: {
				ss << "'" << static_cast<StringParameter*>(it->second)->get() << "'";
				break;
			}
			case INTEGER_PARAMETER: {
				ss << "'" << static_cast<IntegerParameter*>(it->second)->get() << "'";
				break;
			}
			case DOUBLE_PARAMETER: {
				ss << "'" << static_cast<DoubleParameter*>(it->second)->get() << "'";
				break;
			}
			case BOOLEAN_PARAMETER: {
				if (static_cast<BooleanParameter*>(it->second)->get()) {
					ss << "true";
				} else {
					ss << "false";
				}
				break;
			}
		}
		ESINFO(ALWAYS) << ss.str();
	}
}

int Configuration::_getValue(const std::string &parameter, int defaultValue) const
{
	if (_parameters.find(parameter)->second->isSet()) {
		return static_cast<IntegerParameter*>(_parameters.find(parameter)->second)->get();
	} else {
		return defaultValue;
	}
}

long Configuration::_getValue(const std::string &parameter, long defaultValue) const
{
	if (_parameters.find(parameter)->second->isSet()) {
		return static_cast<IntegerParameter*>(_parameters.find(parameter)->second)->get();
	} else {
		return defaultValue;
	}
}

double Configuration::_getValue(const std::string &parameter, double defaultValue) const
{
	if (_parameters.find(parameter)->second->isSet()) {
		return static_cast<DoubleParameter*>(_parameters.find(parameter)->second)->get();
	} else {
		return defaultValue;
	}
}

std::string Configuration::_getValue(const std::string &parameter, const std::string &defaultValue) const
{
	if (_parameters.find(parameter)->second->isSet()) {
		return static_cast<StringParameter*>(_parameters.find(parameter)->second)->get();
	} else {
		return defaultValue;
	}
}

bool Configuration::_getValue(const std::string &parameter, bool defaultValue) const
{
	if (_parameters.find(parameter)->second->isSet()) {
		return static_cast<BooleanParameter*>(_parameters.find(parameter)->second)->get();
	} else {
		return defaultValue;
	}
}

Configuration::~Configuration()
{
	std::map<std::string, Parameter*>::iterator it;
	for (it = _parameters.begin(); it != _parameters.end(); ++it) {
		delete it->second;
	}
}