
#include "meshgenerator.h"

#include "composition/gridgenerator.h"
#include "composition/gridtowergenerator.h"
#include "composition/spheregenerator.h"

#include "mesh/mesh.h"
#include "basis/logging/logging.h"
#include "config/ecf/input/generator.h"

using namespace espreso;

double MeshGenerator::precision = 1e-4;

void MeshGenerator::generate(const InputGeneratorConfiguration &configuration, Mesh &mesh)
{
	switch (configuration.shape) {
	case INPUT_GENERATOR_SHAPE::GRID:
		GridGenerator::generate(configuration.grid, mesh);
		mesh.uniformDecomposition = BlockSettings::uniformDistribution(configuration.grid);
		mesh.preferedDomains = BlockSettings::preferedDomains(configuration.grid);
		break;
	case INPUT_GENERATOR_SHAPE::GRID_TOWER:
		GridTowerGenerator::generate(configuration.grid_tower, mesh);
		mesh.uniformDecomposition = BlockSettings::uniformDistribution(configuration.grid_tower.grids.at(GridTowerGenerator::gridIndex(configuration.grid_tower)));
		mesh.preferedDomains = BlockSettings::preferedDomains(configuration.grid_tower.grids.at(GridTowerGenerator::gridIndex(configuration.grid_tower)));
		break;
	case INPUT_GENERATOR_SHAPE::SPHERE:
		SphereGenerator::generate(configuration.sphere, mesh);
		mesh.uniformDecomposition = BlockSettings::uniformDistribution(configuration.sphere);
		mesh.preferedDomains = BlockSettings::preferedDomains(configuration.sphere);
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Not implemented mesh generator shape";
	}
}



