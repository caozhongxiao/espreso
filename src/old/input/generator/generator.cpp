
#include "generator.h"

#include "../../../config/ecf/input/generator.h"
#include "../../../basis/logging/logging.h"
#include "composition/gridgenerator.h"
#include "composition/gridtowergenerator.h"
#include "composition/spheregenerator.h"

using namespace espreso::input;

double Generator::precision = 1e-4;

void Generator::generate(const InputGeneratorConfiguration &configuration, Mesh &mesh, size_t index, size_t size)
{
	switch (configuration.shape) {
	case INPUT_GENERATOR_SHAPE::GRID:
		GridGenerator::load(configuration.grid, mesh, index, size);
		break;
	case INPUT_GENERATOR_SHAPE::GRID_TOWER:
		GridTowerGenerator::load(configuration.grid_tower, mesh, index, size);
		break;
	case INPUT_GENERATOR_SHAPE::SPHERE:
		SphereGenerator::load(configuration.sphere, mesh, index, size);
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Not implemented generator";
	}
}



