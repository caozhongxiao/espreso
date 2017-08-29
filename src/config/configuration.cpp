
#include "configuration.h"

#include "../basis/logging/logging.h"
#include "../basis/utilities/parser.h"
#include "../mesh/settings/evaluator.h"

#include <iostream>
#include <algorithm>

using namespace espreso;

ECFExpression::ECFExpression()
: evaluator(NULL)
{

}

ECFExpression::~ECFExpression()
{
	if (evaluator) {
		delete evaluator;
	}
}

double ECFExpression::evaluate(const Point &p, double time, double temperature, double pressure, double velocity) const
{
	if (evaluator) {
		return evaluator->evaluate(p, time, temperature, pressure, velocity);
	}
	ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: not created ECFExpression evaluator.";
	return 0;
}

void ECFMetaData::checkdescription(const std::string &name, size_t size) const
{
	if (description.size() != size) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: " << name << " has incorrect number of descriptions.";
	}
}

void ECFMetaData::checkdatatype(const std::string &name, size_t size) const
{
	if (datatype.size() != size) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: " << name << " has incorrect number of datatypes.";
	}
}

void ECFMetaData::checkdefaultvalue(const std::string &name, size_t size) const
{
	if (pattern.size() != size) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: " << name << " has incorrect number of pattern values.";
	}
}

template <typename TType>
static std::vector<TType> getsuffix(size_t start, const std::vector<TType> &data)
{
	if (start < data.size()) {
		return std::vector<TType>(data.begin() + start, data.end());
	}
	return std::vector<TType>();
}

ECFMetaData ECFMetaData::suffix(size_t start) const
{
	return ECFMetaData()
			.setdescription(getsuffix(start, description))
			.setdatatype(getsuffix(start, datatype))
			.setpattern(getsuffix(start, pattern));
}

std::string ECFObject::getValue() const
{
	ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: calling 'get' on ECFObject.";
	return "";
}

bool ECFObject::setValue(const std::string &value)
{
	ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: calling 'set' on ECFObject.";
	return false;
}

ECFParameter* ECFObject::getParameter(const std::string &name)
{
	for (size_t i = 0; i < parameters.size(); i++) {
		if (StringCompare::caseInsensitiveEq(name, parameters[i]->name)) {
			return parameters[i];
		}
	}
	return NULL;
}

ECFParameter* ECFObject::getWithError(const std::string &name)
{
	if (getParameter(name) == NULL) {
		ESINFO(GLOBAL_ERROR) << "ECF ERROR: Object " << this->name << " has no parameter '" << name << "'";
	}
	return getParameter(name);
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

void ECFObject::addSeparator()
{
	registerParameter("SEPARATOR", new ECFSeparator(), ECFMetaData().setdatatype({ ECFDataType::SEPARATOR }));
}

void ECFObject::addSpace()
{
	registerParameter("SPACE", new ECFSeparator(), ECFMetaData().setdatatype({ ECFDataType::SPACE }));
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
	return parameters.back();
}

ECFParameter* ECFObject::registerPatternParameter(ECFParameter *parameter) const
{
	registeredParameters.push_back(parameter);
	parameter->name = metadata.pattern.front();
	parameter->metadata = metadata.suffix(1);
	return parameter;
}

void ECFObject::forEachParameters(std::function<void(const ECFParameter*)> fnc, bool onlyAllowed) const
{
	for (size_t i = 0; i < parameters.size(); i++) {
		if (onlyAllowed && !parameters[i]->metadata.isallowed()) {
			continue;
		}
		if (parameters[i]->isObject()) {
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

