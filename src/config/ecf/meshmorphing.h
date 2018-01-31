
#ifndef SRC_CONFIGURATION_MESHMORPHING_H_
#define SRC_CONFIGURATION_MESHMORPHING_H_

#include "../configuration.h"
#include "../expression.h"
#include "material/coordinatesystem.h"

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
	FIXED,
	OFFSET,
	SCALING,
	TRANSLATION,
	ROTATION
};

struct RBFTargetTransformationConfiguration: public ECFObject {

	MORPHING_TRANSFORMATION transformation;

	ECFExpression offset;
	ECFExpressionVector scaling, translation;
	CoordinateSystemConfiguration coordinate_system;

	bool override;

	RBFTargetTransformationConfiguration(ECFConfiguration *ECFRoot);

protected:
	ECFConfiguration *_ECFRoot;
};

struct RBFTargetConfiguration: public ECFObject {

	MORPHING_RBF_SOLVER solver;

	ECFExpression function;
	double solver_precision;

	std::string target;
	std::map<std::string, RBFTargetTransformationConfiguration> morphers;

	RBFTargetConfiguration(ECFConfiguration *ECFRoot);
};

struct MeshMorphing: public ECFObject {

	MORPHING_TYPE type;
	std::map<std::string, RBFTargetConfiguration> rbf;

	MeshMorphing(ECFConfiguration *ECFRoot);
};

}

#endif /* SRC_CONFIGURATION_MESHMORPHING_H_ */
