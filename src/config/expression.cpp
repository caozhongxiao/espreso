
#include "expression.h"
#include "configuration.hpp"

#include "../basis/utilities/parser.h"
#include "../basis/utilities/utils.h"
#include "../basis/logging/logging.h"
#include "../basis/expression/expression.h"
#include "../basis/evaluator/constevaluator.h"
#include "../basis/evaluator/expressionevaluator.h"
#include "../basis/evaluator/tableinterpolationevaluator.h"

#include "../mesh/store/elementsregionstore.h"
#include "../mesh/store/boundaryregionstore.h"

using namespace espreso;

ECFExpression::ECFExpression(const std::vector<std::string> &variables)
: variables(variables), evaluator(NULL)
{

}

ECFExpression::ECFExpression(const std::vector<std::string> &variables, const std::string &initialValue)
: value(initialValue), variables(variables), evaluator(NULL)
{
	createEvaluator();
}

ECFExpression::ECFExpression(RegionMapBase::RegionIntersection intersection, std::vector<std::string> &regions, const std::map<std::string, ECFExpression> &values)
: variables(values.begin()->second.variables), evaluator(NULL)
{
	std::vector<std::string> filled;
	for (size_t i = 0; i < regions.size(); i++) {
		if (values.at(regions[i]).evaluator != NULL) {
			filled.push_back(regions[i]);
		}
	}
	if (filled.size() == 0) {
		return;
	}

	std::string expression = "(";
	for (size_t i = 0; i + 1 < filled.size(); i++) {
		expression += values.at(filled[i]).evaluator->getEXPRTKForm() + " + ";
	}
	expression += values.at(filled.back()).evaluator->getEXPRTKForm() + ")";

	switch (intersection) {
	case RegionMapBase::RegionIntersection::FIRST:
		value = values.at(filled.front()).value;
		evaluator = values.at(filled.front()).evaluator->copy();
		break;
	case RegionMapBase::RegionIntersection::LAST:
		value = values.at(filled.back()).value;
		evaluator = values.at(filled.back()).evaluator->copy();
		break;
	case RegionMapBase::RegionIntersection::AVERAGE:
		value = expression + " / " + std::to_string(filled.size());
		createEvaluator();
		break;
	case RegionMapBase::RegionIntersection::SUM:
		value = expression;
		createEvaluator();
		break;
	case RegionMapBase::RegionIntersection::ERROR:
		ESINFO(ERROR) << "The following regions intersect: " << regions;
		break;
	default:
		ESINFO(ERROR) << "ESPRESO internal error: not implemented region intersection type.";
	}
}

ECFExpression::~ECFExpression()
{
	if (evaluator) {
		delete evaluator;
	}
}

ECFExpression::ECFExpression(const ECFExpression &other)
{
	value = other.value;
	variables = other.variables;
	evaluator = NULL;
	if (other.evaluator != NULL) {
		evaluator = other.evaluator->copy();
	}
}

ECFExpression& ECFExpression::operator=(const ECFExpression &other)
{
	if (this != &other) {
		value = other.value;
		variables = other.variables;
		if (evaluator != NULL) {
			delete evaluator;
			evaluator = NULL;
		}
		if (other.evaluator != NULL) {
			evaluator = other.evaluator->copy();
		}
	}
	return *this;
}

bool ECFExpression::createEvaluator()
{
	if (evaluator != NULL) {
		delete evaluator;
	}
	if (StringCompare::contains(this->value, { "TABULAR" })) {
		std::string value = Parser::strip(this->value.substr(this->value.find_first_of("[")));
		value = value.substr(1, value.size() - 3);
		std::vector<std::string> lines = Parser::split(value, ";");
		std::vector<std::pair<double, double> > table;

		for (size_t i = 0; i < lines.size(); i++) {
			if (lines[i].size() == 0) {
				continue;
			}
			std::vector<std::string> line = Parser::split(lines[i], ",");
			if (line.size() != 2) {
				ESINFO(GLOBAL_ERROR) << "Invalid TABULAR data: " << value;
			}
			table.push_back(std::make_pair(std::stod(line[0]), std::stod(line[1])));
		}
		evaluator = new TableInterpolationEvaluator(table);
		return true;
	}
	if (StringCompare::contains(this->value, ExpressionEvaluator::variables())) {
		if (Expression::isValid(this->value, variables)) {
			evaluator = new ExpressionEvaluator(this->value);
			return true;
		}
	} else {
		if (Expression::isValid(this->value, variables)) {
			Expression expr(this->value, {});
			evaluator = new ConstEvaluator(expr.evaluate());
			return true;
		}
	}

	return false;
}


void ECFExpressionVector::init()
{
	REGISTER(x, ECFMetaData()
			.setdescription({ "x-direction." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.allowonly([&] () { return dimension != DIMENSION::Z; }));
	REGISTER(y, ECFMetaData()
			.setdescription({ "y-direction." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.allowonly([&] () { return dimension == DIMENSION::D2 || dimension == DIMENSION::D3; }));
	REGISTER(z, ECFMetaData()
			.setdescription({ "z-direction." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.allowonly([&] () { return dimension == DIMENSION::Z || dimension == DIMENSION::D3; }));
}

ECFExpressionVector::ECFExpressionVector(DIMENSION dimension, const std::vector<std::string> &variables)
: x(variables), y(variables), z(variables), dimension(dimension)
{
	init();
}

ECFExpressionVector::ECFExpressionVector(DIMENSION dimension, const std::vector<std::string> &variables, const std::string &initialValue)
: x(variables, initialValue), y(variables, initialValue), z(variables, initialValue), dimension(dimension)
{
	init();
}

ECFExpressionVector::ECFExpressionVector(RegionMapBase::RegionIntersection intersection, std::vector<std::string> &regions, const std::map<std::string, ECFExpressionVector> &values)
: x(intersection, regions, getComponent(values, [] (std::map<std::string, ECFExpressionVector>::const_iterator it) { return it->second.x; })),
  y(intersection, regions, getComponent(values, [] (std::map<std::string, ECFExpressionVector>::const_iterator it) { return it->second.y; })),
  z(intersection, regions, getComponent(values, [] (std::map<std::string, ECFExpressionVector>::const_iterator it) { return it->second.z; })),
  dimension(values.begin()->second.dimension)
{

}

ECFExpressionVector::ECFExpressionVector(RegionMapBase::RegionIntersection intersection, std::vector<std::string> &regions, const std::map<std::string, ECFExpressionOptionalVector> &values)
: x(intersection, regions, getComponent(values, [] (std::map<std::string, ECFExpressionOptionalVector>::const_iterator it) { return it->second.x; })),
  y(intersection, regions, getComponent(values, [] (std::map<std::string, ECFExpressionOptionalVector>::const_iterator it) { return it->second.y; })),
  z(intersection, regions, getComponent(values, [] (std::map<std::string, ECFExpressionOptionalVector>::const_iterator it) { return it->second.z; })),
  dimension(values.begin()->second.dimension)
{

}

ECFExpressionOptionalVector::ECFExpressionOptionalVector(DIMENSION dimension, const std::vector<std::string> &variables)
: ECFExpressionVector(dimension, variables), all(variables)
{
	REGISTER(all, ECFMetaData()
			.setdescription({ "all-directions." })
			.setdatatype({ ECFDataType::EXPRESSION }));
}

ECFExpressionOptionalVector::ECFExpressionOptionalVector(RegionMapBase::RegionIntersection intersection, std::vector<std::string> &regions, const std::map<std::string, ECFExpressionOptionalVector> &values)
: ECFExpressionVector(intersection, regions, values),
  all(intersection, regions, getComponent(values, [] (std::map<std::string, ECFExpressionOptionalVector>::const_iterator it) { return it->second.all; }))
{

}


