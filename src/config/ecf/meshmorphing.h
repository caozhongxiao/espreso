
#ifndef SRC_CONFIGURATION_MESHMORPHING_H_
#define SRC_CONFIGURATION_MESHMORPHING_H_

#include "../configuration.h"
#include "../expression.h"

namespace espreso {

class ECFConfiguration;

enum class MORPHING_TYPE {
	NONE = 0,
	RBF = 1
};

enum class MORPHING_RBF_SOLVER {
	ITERATIVE = 0,
	DIRECT = 1
};

enum class MORPHING_TRANSFORMATION {
	FIXED = 0,
	OFFSET = 1,
	ROTATION = 2,
	TRANSLATION = 3
};

struct RBFTargetTransformation : public ECFObject {

	MORPHING_TRANSFORMATION transformation;

	ECFExpression offset;
	ECFExpressionVector translation;

	bool overriding;

	RBFTargetTransformation(ECFConfiguration *ECFRoot);

protected:
	ECFConfiguration *_ECFRoot;
};

struct RBFTarget : public ECFObject {

	MORPHING_RBF_SOLVER solver;

	ECFExpression function;
	double solver_precision;

	std::map<std::string, RBFTargetTransformation> targets;

	RBFTarget(ECFConfiguration *ECFRoot);
};

struct MeshMorphing: public ECFObject {

	MORPHING_TYPE type;
	std::map<std::string, RBFTarget> rbf;

	MeshMorphing(ECFConfiguration *ECFRoot);
};

}

#endif /* SRC_CONFIGURATION_MESHMORPHING_H_ */
