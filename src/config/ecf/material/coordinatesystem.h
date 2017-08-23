
#ifndef SRC_CONFIG_ECF_MATERIAL_COORDINATESYSTEM_H_
#define SRC_CONFIG_ECF_MATERIAL_COORDINATESYSTEM_H_

#include "../../configuration.h"

namespace espreso {

struct CoordinateSystemConfiguration: public ECFObject {

	enum class TYPE {
		CARTESIAN,
		CYLINDRICAL,
		SPHERICAL
	};

	TYPE type;

	ECFExpression rotation_x, rotation_y, rotation_z;
	ECFExpression center_x, center_y, center_z;

	CoordinateSystemConfiguration();
};


}

#endif /* SRC_CONFIG_ECF_MATERIAL_COORDINATESYSTEM_H_ */
