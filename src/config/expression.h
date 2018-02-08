
#ifndef SRC_CONFIG_EXPRESSION_H_
#define SRC_CONFIG_EXPRESSION_H_

#include "configuration.h"
#include <functional>

namespace espreso {

class Evaluator;
struct BoundaryRegionStore;
struct ElementsRegionStore;

struct ECFExpression {
	std::string value;
	std::vector<std::string> variables;
	Evaluator *evaluator;

	ECFExpression(const std::vector<std::string> &variables);
	ECFExpression(const std::vector<std::string> &variables, const std::string &initialValue);
	ECFExpression(const ECFExpression &other);
	ECFExpression& operator=(const ECFExpression &other);
	~ECFExpression();

	ECFExpression(std::vector<std::string> &regions, const std::map<std::string, ECFExpression> &values);

	bool createEvaluator();
};

struct ECFExpressionOptionalVector;

struct ECFExpressionVector: public ECFObject {
	ECFExpression x, y, z;
	DIMENSION dimension;

	ECFExpressionVector(DIMENSION dimension, const std::vector<std::string> &variables);
	ECFExpressionVector(DIMENSION dimension, const std::vector<std::string> &variables, const std::string &initialValue);

	ECFExpressionVector(std::vector<std::string> &regions, const std::map<std::string, ECFExpressionVector> &values);
	ECFExpressionVector(std::vector<std::string> &regions, const std::map<std::string, ECFExpressionOptionalVector> &values);

protected:
	template <typename TValue>
	std::map<std::string, ECFExpression> getComponent(
			const std::map<std::string, TValue> &values,
			std::function<ECFExpression(typename std::map<std::string, TValue>::const_iterator)> get)
	{
		std::map<std::string, ECFExpression> component;
		for (auto it = values.begin(); it != values.end(); ++it) {
			component.emplace(std::piecewise_construct, std::forward_as_tuple(it->first), std::forward_as_tuple(get(it)));
		}
		return component;
	}

	void init();
};

struct ECFExpressionOptionalVector: public ECFExpressionVector {
	ECFExpression all;

	ECFExpressionOptionalVector(DIMENSION dimension, const std::vector<std::string> &variables);

	ECFExpressionOptionalVector(std::vector<std::string> &regions, const std::map<std::string, ECFExpressionOptionalVector> &values);
};

}

#endif /* SRC_CONFIG_EXPRESSION_H_ */
