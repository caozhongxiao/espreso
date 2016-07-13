
#ifndef INPUT_MESHGENERATOR_CONFIGURATION_PARAMETER_H_
#define INPUT_MESHGENERATOR_CONFIGURATION_PARAMETER_H_

#include <string>
#include <stdlib.h>
#include <iostream>
#include <sstream>

#include "../logging/logging.h"
#include "../options/options.h"

namespace espreso {
namespace input {

enum DataType {
	STRING_PARAMETER,
	INTEGER_PARAMETER,
	LONG_PARAMETER,
	SIZE_PARAMETER,
	DOUBLE_PARAMETER,
	BOOLEAN_PARAMETER
};

enum HelpFlags {
	INGNORE_IN_HELP,
	WRITE_TO_HELP
};


struct Description {
	DataType type;
	std::string name;
	void* value;
	std::string description;
	std::vector<std::string> options;
	bool writeToHelp;

	Description(std::string name, int &defaultValue, std::string description, HelpFlags writeToHelp = INGNORE_IN_HELP)
	: type(INTEGER_PARAMETER), name(name), value(&defaultValue), description(description), writeToHelp(writeToHelp) { };

	Description(std::string name, long &defaultValue, std::string description, HelpFlags writeToHelp = INGNORE_IN_HELP)
	: type(INTEGER_PARAMETER), name(name), value(&defaultValue), description(description), writeToHelp(writeToHelp) { };

	Description(std::string name, size_t &defaultValue, std::string description, HelpFlags writeToHelp = INGNORE_IN_HELP)
	: type(INTEGER_PARAMETER), name(name), value(&defaultValue), description(description), writeToHelp(writeToHelp) { };

	Description(std::string name, double &defaultValue, std::string description, HelpFlags writeToHelp = INGNORE_IN_HELP)
	: type(DOUBLE_PARAMETER), name(name), value(&defaultValue), description(description), writeToHelp(writeToHelp) { };

	Description(std::string name, std::string &defaultValue, std::string description, HelpFlags writeToHelp = INGNORE_IN_HELP)
	: type(STRING_PARAMETER), name(name), value(&defaultValue), description(description), writeToHelp(writeToHelp) { };

	Description(std::string name, bool &defaultValue, std::string description, HelpFlags writeToHelp = INGNORE_IN_HELP)
	: type(BOOLEAN_PARAMETER), name(name), value(&defaultValue), description(description), writeToHelp(writeToHelp) { };


	Description(std::string name, int &defaultValue, std::string description, std::vector<std::string> options, HelpFlags writeToHelp = INGNORE_IN_HELP)
	: type(INTEGER_PARAMETER), name(name), value(&defaultValue), description(description), options(options), writeToHelp(writeToHelp) { };

	Description(std::string name, long &defaultValue, std::string description, std::vector<std::string> options, HelpFlags writeToHelp = INGNORE_IN_HELP)
	: type(INTEGER_PARAMETER), name(name), value(&defaultValue), description(description), options(options), writeToHelp(writeToHelp) { };

	Description(std::string name, size_t &defaultValue, std::string description, std::vector<std::string> options, HelpFlags writeToHelp = INGNORE_IN_HELP)
	: type(INTEGER_PARAMETER), name(name), value(&defaultValue), description(description), options(options), writeToHelp(writeToHelp) { };

	Description(std::string name, double &defaultValue, std::string description, std::vector<std::string> options, HelpFlags writeToHelp = INGNORE_IN_HELP)
	: type(DOUBLE_PARAMETER), name(name), value(&defaultValue), description(description), options(options), writeToHelp(writeToHelp) { };

	Description(std::string name, std::string &defaultValue, std::string description, std::vector<std::string> options, HelpFlags writeToHelp = INGNORE_IN_HELP)
	: type(STRING_PARAMETER), name(name), value(&defaultValue), description(description), options(options), writeToHelp(writeToHelp) { };

	Description(std::string name, bool &defaultValue, std::string description, std::vector<std::string> options, HelpFlags writeToHelp = INGNORE_IN_HELP)
	: type(BOOLEAN_PARAMETER), name(name), value(&defaultValue), description(description), options(options), writeToHelp(writeToHelp) { };
};

struct CaseInsensitiveCompare {

	bool operator()(const std::string &p1, const std::string &p2) const
	{
		return std::lexicographical_compare(p1.begin(), p1.end(), p2.begin(), p2.end(), caseInsensitive);
	}

	static bool caseInsensitive(const char &c1, const char &c2)
	{
		return std::tolower(c1) < std::tolower(c2);
	}

	static bool equals(const char &c1, const char &c2)
	{
		return std::tolower(c1) == std::tolower(c2);
	}
};

class Parameter {

public:
	bool match(const std::string &line) const
	{
		std::string param = line.substr(0, line.find_first_of(" ="));
		return param.size() == _name.size() && std::equal(param.begin(), param.end(), _name.begin(), CaseInsensitiveCompare::equals);
	}

	const std::string& name() const
	{
		return _name;
	}

	const std::string& description() const
	{
		return _description;
	}

	DataType type() const
	{
		return _type;
	}

	bool isSet() const
	{
		return _set;
	}

	void reset(bool value)
	{
		_set = value;
	}

	virtual void set(const std::string &line) =0;
	virtual Parameter* copy() =0;

	virtual ~Parameter() {};

protected:
	Parameter(DataType type, std::string name, std::string description)
		: _type(type), _delimiter("="), _name(name), _description(description), _set(false) {};

	std::string value(const std::string &line)
	{
		size_t pos = line.find(_delimiter);
		if (pos == std::string::npos) {
			ESINFO(ERROR) << "Incorrect format of " << _name << ". Use " << _name << _delimiter << "value.";
		}
		std::string val = line.substr(pos + 1);
		val.erase(0, val.find_first_not_of(" "));
		if (val[val.size() - 1] == ' ' &&  val.find_last_of(" ") != std::string::npos) {
			return val.erase(val.find_last_of(" "));
		} else {
			return val;
		}
	}

	std::string _delimiter;
	std::string _name;
	std::string _description;
	DataType _type;
	bool _set;

};

class StringParameter : public Parameter {

public:
	StringParameter(std::string name, std::string description)
		:Parameter(STRING_PARAMETER, name, description), _value("") {};

	void set(const std::string &line)
	{
		_value = value(line);
		if (!_value.size()) {
			ESINFO(ERROR) << "Empty parameter " << _name << ".";
		}
		_set = true;
	}

	const std::string& get() const { return _value; };

	Parameter* copy()
	{
		return new StringParameter(*this);
	}

private:
	std::string _value;
};

class IntegerParameter : public Parameter {

public:
	IntegerParameter(std::string name, std::string description)
		:Parameter(INTEGER_PARAMETER, name, description), _value(0) {};

	void set(const std::string &line)
	{
		std::stringstream ss(value(line));
		ss >> _value;
		_set = true;
	}

	const eslocal& get() const { return _value; };

	Parameter* copy()
	{
		return new IntegerParameter(*this);
	}

private:
	eslocal _value;
};

class DoubleParameter : public Parameter {

public:
	DoubleParameter(std::string name, std::string description)
		:Parameter(DOUBLE_PARAMETER, name, description), _value(0) {};

	void set(const std::string &line)
	{
		std::stringstream ss(value(line));
		ss >> _value;
		_set = true;
	}

	const double& get() const { return _value; };

	Parameter* copy()
	{
		return new DoubleParameter(*this);
	}

private:
	double _value;
};

class BooleanParameter : public Parameter {

public:
	BooleanParameter(std::string name, std::string description)
		:Parameter(BOOLEAN_PARAMETER, name, description), _value(false) {};

	void set(const std::string &line)
	{
		size_t pos = line.find(_delimiter);
		if (pos == std::string::npos) {
			_value = true;
		} else {
			if (value(line).compare("0") == 0) {
				_value = false;
			} else {
				_value = true;
			}
		}
		_set = true;
	}

	const bool& get() const { return _value; };

	Parameter* copy()
	{
		return new BooleanParameter(*this);
	}

private:
	bool _value;
};

}
}


#endif /* INPUT_MESHGENERATOR_CONFIGURATION_PARAMETER_H_ */