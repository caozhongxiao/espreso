
#include "expression.h"
#include "configuration.hpp"

#include "../basis/utilities/parser.h"
#include "../basis/logging/logging.h"
#include "../basis/expression/expression.h"
#include "../basis/evaluator/constevaluator.h"
#include "../basis/evaluator/expressionevaluator.h"
#include "../basis/evaluator/tableinterpolationevaluator.h"

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

ECFExpression::ECFExpression(const ECFExpression &other)
{
	value = other.value;
	if (other.evaluator != NULL) {
		evaluator = other.evaluator->copy();
	}
}

ECFExpression& ECFExpression::operator=(const ECFExpression &other)
{
	if (this != &other) {
		value = other.value;
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

bool ECFExpression::createEvaluator(const std::vector<std::string> &variables)
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

espreso::ECFExpressionVector::ECFExpressionVector(DIMENSION dimension, bool fillWithZeros)
: dimension(dimension)
{
	REGISTER(x, ECFMetaData()
			.setdescription({ "x-direction." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setboundaryconditionvariables()
			.allowonly([&] () { return dimension != DIMENSION::Z; }));
	REGISTER(y, ECFMetaData()
			.setdescription({ "y-direction." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setboundaryconditionvariables()
			.allowonly([&] () { return dimension == DIMENSION::D2 || dimension == DIMENSION::D3; }));
	REGISTER(z, ECFMetaData()
			.setdescription({ "z-direction." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setboundaryconditionvariables()
			.allowonly([&] () { return dimension == DIMENSION::Z || dimension == DIMENSION::D3; }));

	if (fillWithZeros) {
		x.value = "0";
		y.value = "0";
		z.value = "0";
		x.createEvaluator(ECFMetaData().setboundaryconditionvariables().variables);
		y.createEvaluator(ECFMetaData().setboundaryconditionvariables().variables);
		z.createEvaluator(ECFMetaData().setboundaryconditionvariables().variables);
	}
}

espreso::ECFExpressionOptionalVector::ECFExpressionOptionalVector(DIMENSION dimension)
: ECFExpressionVector(dimension, false)
{
	REGISTER(all, ECFMetaData()
			.setdescription({ "all-directions." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setboundaryconditionvariables());
}


