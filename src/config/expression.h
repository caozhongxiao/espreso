
#ifndef SRC_CONFIG_EXPRESSION_H_
#define SRC_CONFIG_EXPRESSION_H_

#include "config/configuration.h"
#include <functional>

namespace espreso {

class Evaluator;
struct BoundaryRegionStore;
struct ElementsRegionStore;

struct ECFExpression {
	std::string value;
	std::vector<std::string> variables;
	Evaluator *evaluator;

	bool isSet() const { return value.size(); }

	ECFExpression(const std::vector<std::string> &variables);
	ECFExpression(const std::vector<std::string> &variables, const std::string &initialValue);
	ECFExpression(const ECFExpression &other);
	ECFExpression& operator=(const ECFExpression &other);
	~ECFExpression();

	ECFExpression(RegionMapBase::RegionIntersection intersection, std::vector<std::string> &regions, const std::map<std::string, ECFExpression> &values);

	bool createEvaluator();
};

struct ECFExpressionOptionalVector;

struct ECFExpressionVector: public ECFObject {
	ECFExpression data[3];
	ECFExpression &x = data[0], &y = data[1], &z = data[2];
	DIMENSION dimension;

	ECFExpressionVector(const ECFExpressionVector &other);
	ECFExpressionVector& operator=(const ECFExpressionVector &other);
	ECFExpressionVector(DIMENSION dimension, const std::vector<std::string> &variables);
	ECFExpressionVector(DIMENSION dimension, const std::vector<std::string> &variables, const std::string &initialValue);

	ECFExpressionVector(RegionMapBase::RegionIntersection intersection, std::vector<std::string> &regions, const std::map<std::string, ECFExpressionVector> &values);
	ECFExpressionVector(RegionMapBase::RegionIntersection intersection, std::vector<std::string> &regions, const std::map<std::string, ECFExpressionOptionalVector> &values);

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

	ECFExpressionOptionalVector(RegionMapBase::RegionIntersection intersection, std::vector<std::string> &regions, const std::map<std::string, ECFExpressionOptionalVector> &values);
};

}

#endif /* SRC_CONFIG_EXPRESSION_H_ */
