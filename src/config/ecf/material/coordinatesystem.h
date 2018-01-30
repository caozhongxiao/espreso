
#ifndef SRC_CONFIG_ECF_MATERIAL_COORDINATESYSTEM_H_
#define SRC_CONFIG_ECF_MATERIAL_COORDINATESYSTEM_H_

#include "../../configuration.h"
#include "../../expression.h"

namespace espreso {

struct CoordinateSystemConfiguration: public ECFObject {

	enum class TYPE {
		CARTESIAN,
		CYLINDRICAL,
		SPHERICAL
	};

	TYPE type;
	DIMENSION dimension;

	ECFExpressionVector rotation;
	ECFExpressionVector center;

	CoordinateSystemConfiguration();
};


}

#endif /* SRC_CONFIG_ECF_MATERIAL_COORDINATESYSTEM_H_ */
