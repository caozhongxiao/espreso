
#ifndef SRC_CONFIG_ECF_INPUT_GENERATOR_H_
#define SRC_CONFIG_ECF_INPUT_GENERATOR_H_

#include "config/configuration.h"

#include "grid.h"
#include "gridtower.h"
#include "sphere.h"

namespace espreso {

enum class INPUT_GENERATOR_SHAPE {
	GRID,
	GRID_TOWER,
	SPHERE
};

struct InputGeneratorConfiguration: public ECFObject {

	INPUT_GENERATOR_SHAPE shape;

	GridGeneratorConfiguration grid;
	GridTowerGeneratorConfiguration grid_tower;
	SphereGeneratorConfiguration sphere;

	InputGeneratorConfiguration();
};

}



#endif /* SRC_CONFIG_ECF_INPUT_GENERATOR_H_ */
