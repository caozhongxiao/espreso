
#ifndef SRC_CONFIGURATION_MESHMORPHING_H_
#define SRC_CONFIGURATION_MESHMORPHING_H_

#include <config/holders/expression.h>
#include "config/configuration.h"
#include "config/ecf/material/coordinatesystem.h"

namespace espreso {

class ECFRoot;

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

	RBFTargetTransformationConfiguration(ECFRoot *ECFRoot);
protected:
	ECFRoot *_ECFRoot;
};

struct ExternalFFDConfiguration: public ECFObject {

	std::string path;
	std::map<std::string, RBFTargetTransformationConfiguration> morphers;

	ExternalFFDConfiguration(ECFRoot *ECFRoot);
};

struct RBFTargetConfiguration: public ECFObject {

	MORPHING_RBF_SOLVER solver;

	ECFExpression function;
	double solver_precision;
	int solver_max_iter;

	std::string target;
	std::map<std::string, RBFTargetTransformationConfiguration> morphers;

	ExternalFFDConfiguration external_ffd;

	RBFTargetConfiguration(ECFRoot *ECFRoot);
};

struct MeshMorphing: public ECFObject {

	MORPHING_TYPE type;
	std::map<std::string, RBFTargetConfiguration> rbf;

	MeshMorphing(ECFRoot *ECFRoot);
};

}

#endif /* SRC_CONFIGURATION_MESHMORPHING_H_ */
