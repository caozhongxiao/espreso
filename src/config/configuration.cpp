
#include "configuration.h"
#include "conditions.h"

#include "basis/logging/logging.h"
#include "basis/utilities/parser.h"
#include <iostream>
#include <algorithm>

using namespace espreso;

void ECFParameter::defaultName()
{
	if (!metadata.name.size()) {
		metadata.name = name;

		for (size_t i = 0, up = 1; i < metadata.name.size(); i++) {
			if (up && 'a' <= metadata.name[i] && metadata.name[i] <= 'z') {
				metadata.name[i] += 'A' - 'a';
			}
			up = 0;
			if (metadata.name[i] == '_') {
				metadata.name[i] = ' ';
				up = 1;
			}
		}
	}
}

bool ECFParameter::isvisible()
{
	return metadata.condition->evaluate();
}


bool ECFParameter::setValue(const std::string &value)
{
	if (_setValue(value)) {
		for (size_t i = 0; i < _setValueListeners.size(); i++) {
			_setValueListeners[i](value);
		}
		return true;
	}
	return false;
}

ECFParameter* ECFParameter::_triggerParameterGet(ECFParameter* parameter)
{
	if (parameter != NULL) {
		for (size_t i = 0; i < _parameterGetListeners.size(); i++) {
			_parameterGetListeners[i](parameter->name);
		}
	}
	return parameter;
}

ECFParameter* ECFParameter::getParameter(const std::string &name)
{
	return _triggerParameterGet(_getParameter(name));
}

ECFParameter* ECFParameter::getParameter(const char* name)
{
	return _triggerParameterGet(_getParameter(std::string(name)));
}

ECFParameter* ECFParameter::getParameter(const void* data)
{
	return _triggerParameterGet(_getParameter(data));
}

void ECFParameter::addListener(Event event, std::function<void(const std::string &value)> listener)
{
	switch (event) {
	case Event::VALUE_SET:
		_setValueListeners.push_back(listener);
		break;
	case Event::PARAMETER_GET:
		_parameterGetListeners.push_back(listener);
	}
}

ECFParameter* ECFParameter::registerAdditionalParameter(ECFParameter* parameter)
{
	ESINFO(ERROR) << "ESPRESO internal error: parameter '" << name << "' cannot have additional parameters.";
	return NULL;
}

std::string ECFObject::getValue() const
{
	ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: calling 'get' on ECFObject.";
	return "";
}

bool ECFObject::_setValue(const std::string &value)
{
	ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: calling 'set' on ECFObject.";
	return false;
}

ECFParameter* ECFObject::_getParameter(const std::string &name)
{
	for (size_t i = 0; i < parameters.size(); i++) {
		if (StringCompare::caseInsensitiveEq(name, parameters[i]->name)) {
			return parameters[i];
		}
	}
	return NULL;
}

ECFParameter* ECFObject::_getParameter(const void* data)
{
	for (size_t i = 0; i < parameters.size(); i++) {
		if (data == parameters[i]->data()) {
			return parameters[i];
		}
	}
	return NULL;
}

ECFParameter* ECFObject::getWithError(const std::string &name)
{
	if (_getParameter(name) == NULL) {
		ESINFO(GLOBAL_ERROR) << "ECF ERROR: Object " << this->name << " has no parameter '" << name << "'";
	}
	return _getParameter(name);
}

ECFParameter* ECFObject::registerAdditionalParameter(ECFParameter* parameter)
{
	parameters.push_back(parameter);
	return parameter;
}

void ECFObject::dropParameter(ECFParameter *parameter)
{
	for (size_t i = 0; i < parameters.size(); i++) {
		if (parameters[i] == parameter) {
			parameters.erase(parameters.begin() + i);
			return;
		}
	}
	ESINFO(GLOBAL_ERROR) << "ECF ERROR: Object cannot drop parameter '" << parameter->name << "'";
}

void ECFObject::dropAllParameters()
{
    parameters.clear();
}

ECFParameter* ECFObject::addSeparator()
{
	registerParameter("SEPARATOR", new ECFSeparator(), ECFMetaData().setdatatype({ ECFDataType::SEPARATOR }));
	return parameters.back();
}

ECFParameter* ECFObject::addSpace()
{
	registerParameter("SPACE", new ECFSeparator(), ECFMetaData().setdatatype({ ECFDataType::SPACE }));
	return parameters.back();
}

void ECFObject::moveLastBefore(const std::string &name)
{
	ECFParameter *last = parameters.back();
	auto position = std::find(parameters.begin(), parameters.end(), getWithError(name));
	parameters.pop_back();
	parameters.insert(position, last);
}

ECFParameter* ECFObject::registerParameter(const std::string &name, ECFParameter *parameter, const ECFMetaData &metadata)
{
	registeredParameters.push_back(parameter);
	parameters.push_back(parameter);
	parameters.back()->name = name;
	parameters.back()->metadata = metadata;
	parameters.back()->defaultName();
	return parameters.back();
}

ECFParameter* ECFObject::registerPatternParameter(ECFParameter *parameter) const
{
	registeredParameters.push_back(parameter);
	parameter->name = metadata.pattern.front();
	parameter->metadata = metadata.suffix(1);
	registeredParameters.back()->defaultName();
	return parameter;
}

void ECFObject::forEachParameters(std::function<void(const ECFParameter*)> fnc, bool onlyAllowed) const
{
	for (size_t i = 0; i < parameters.size(); i++) {
		if (onlyAllowed && !parameters[i]->metadata.isallowed()) {
			continue;
		}
		if (parameters[i]->isObject()) {
			fnc(parameters[i]);
			dynamic_cast<const ECFObject*>(parameters[i])->forEachParameters(fnc, onlyAllowed);
		}
		if (parameters[i]->isValue()) {
			fnc(parameters[i]);
		}
	}
}

void ECFObject::forEachParameters(std::function<void(ECFParameter*)> fnc, bool onlyAllowed)
{
	for (size_t i = 0; i < parameters.size(); i++) {
		if (onlyAllowed && !parameters[i]->metadata.isallowed()) {
			continue;
		}
		if (parameters[i]->isObject()) {
			dynamic_cast<ECFObject*>(parameters[i])->forEachParameters(fnc, onlyAllowed);
		}
		if (parameters[i]->isValue()) {
			fnc(parameters[i]);
		}
	}
}

ECFObject::~ECFObject()
{
	for (size_t i = 0; i < registeredParameters.size(); i++) {
		delete registeredParameters[i];
	}
}

