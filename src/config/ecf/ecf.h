
#ifndef SRC_CONFIG_ECF_ECF_H_
#define SRC_CONFIG_ECF_ECF_H_

#include "environment.h"
#include "output.h"
#include "decomposer.h"

#include "input/input.h"
#include "input/generator.h"

#include "physics/physics.h"
#include "physics/advectiondiffusion.h"
#include "physics/structuralmechanics.h"

#include "../reader/reader.h"

namespace espreso {

struct ECFConfiguration: public ECFObject {

	// Environment has to be created first!
	Environment environment;

	std::map<size_t, std::string> default_args;
	std::map<std::string, std::string> variables;

	INPUT_FORMAT input;
	PHYSICS physics;

	InputConfiguration workbench, openfoam, esdata;
	InputGeneratorConfiguration generator;

	AdvectionDiffusion2DConfiguration advection_diffusion_2d;
	AdvectionDiffusion3DConfiguration advection_diffusion_3d;
	StructuralMechanics2DConfiguration structural_mechanics_2d;
	StructuralMechanics3DConfiguration structural_mechanics_3d;

	OutputConfiguration output;

	DecomposerConfiguration decomposer;

	void init();
	ECFConfiguration() { init(); }
	ECFConfiguration(const std::string &file);
	ECFConfiguration(int *argc, char ***argv);
};

}

#endif /* SRC_CONFIG_ECF_ECF_H_ */
