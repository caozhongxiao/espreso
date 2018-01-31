
#include "expression.h"
#include "configuration.hpp"

#include "../basis/utilities/parser.h"
#include "../basis/logging/logging.h"
#include "../basis/expression/expression.h"
#include "../basis/evaluator/constevaluator.h"
#include "../basis/evaluator/expressionevaluator.h"
#include "../basis/evaluator/tableinterpolationevaluator.h"

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
	if (Expression::isValid(this->value, variables)) {
		evaluator = new ExpressionEvaluator(this->value);
		return true;
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

ECFExpressionOptionalVector::ECFExpressionOptionalVector(DIMENSION dimension, const std::vector<std::string> &variables)
: ECFExpressionVector(dimension, variables), all(variables)
{
	REGISTER(all, ECFMetaData()
			.setdescription({ "all-directions." })
			.setdatatype({ ECFDataType::EXPRESSION }));
}


