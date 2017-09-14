
#ifndef SRC_CONFIG_ECF_MATERIAL_COORDINATESYSTEM_H_
#define SRC_CONFIG_ECF_MATERIAL_COORDINATESYSTEM_H_

#include "../../configuration.h"

namespace espreso {

enum class DIMENSION {
	D2,
	D3
};

struct CoordinateSystemConfiguration: public ECFObject {

	enum class TYPE {
		CARTESIAN,
		CYLINDRICAL,
		SPHERICAL
	};

	TYPE type;
	DIMENSION dimension;

	ECFExpression rotation_x, rotation_y, rotation_z;
	ECFExpression center_x, center_y, center_z;

	CoordinateSystemConfiguration();
};


}

#endif /* SRC_CONFIG_ECF_MATERIAL_COORDINATESYSTEM_H_ */
