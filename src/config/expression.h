
#ifndef SRC_CONFIG_EXPRESSION_H_
#define SRC_CONFIG_EXPRESSION_H_

#include "configuration.h"

namespace espreso {

class Evaluator;

struct ECFExpression {
	std::string value;
	std::vector<std::string> variables;
	Evaluator *evaluator;

	ECFExpression(const std::vector<std::string> &variables);
	ECFExpression(const std::vector<std::string> &variables, const std::string &initialValue);
	ECFExpression(const ECFExpression &other);
	ECFExpression& operator=(const ECFExpression &other);
	~ECFExpression();

	bool createEvaluator();
};

struct ECFExpressionVector: public ECFObject {
	ECFExpression x, y, z;
	DIMENSION dimension;

	ECFExpressionVector(DIMENSION dimension, const std::vector<std::string> &variables);
	ECFExpressionVector(DIMENSION dimension, const std::vector<std::string> &variables, const std::string &initialValue);

protected:
	void init();
};

struct ECFExpressionOptionalVector: public ECFExpressionVector {
	ECFExpression all;

	ECFExpressionOptionalVector(DIMENSION dimension, const std::vector<std::string> &variables);
};

}

#endif /* SRC_CONFIG_EXPRESSION_H_ */
