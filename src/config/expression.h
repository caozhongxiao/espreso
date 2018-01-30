
#ifndef SRC_CONFIG_EXPRESSION_H_
#define SRC_CONFIG_EXPRESSION_H_

#include "configuration.h"

namespace espreso {

class Evaluator;

struct ECFExpression {
	std::string value;
	Evaluator *evaluator;

	ECFExpression();
	ECFExpression(const ECFExpression &other);
	ECFExpression& operator=(const ECFExpression &other);
	~ECFExpression();

	bool createEvaluator(const std::vector<std::string> &variables);
};

struct ECFExpressionVector: public ECFObject {
	ECFExpression x, y, z;
	DIMENSION dimension;

	ECFExpressionVector(DIMENSION dimension, bool fillWithZeros);
};

struct ECFExpressionOptionalVector: public ECFExpressionVector {
	ECFExpression all;

	ECFExpressionOptionalVector(DIMENSION dimension);
};

}

#endif /* SRC_CONFIG_EXPRESSION_H_ */
