
#ifndef SRC_CONFIG_ECF_ECF_H_
#define SRC_CONFIG_ECF_ECF_H_

#include "environment.h"
#include "output.h"
#include "decomposer.h"

#include "pythontestgenerator.h"

#include "input/input.h"
#include "input/generator.h"
#include "input/feti4ilibrary.h"

#include "physics/physics.h"
#include "physics/structuralmechanics.h"

#include "../reader/reader.h"
#include "physics/heattransfer.h"

namespace espreso {

struct ECFConfiguration: public ECFObject {

	// Environment has to be created first!
	Environment environment;

	PythonTestGenerator python_test_generator;

	std::map<size_t, std::string> default_args;
	std::map<std::string, std::string> variables;

	INPUT_FORMAT input;
	PHYSICS physics;

	InputConfiguration workbench, openfoam, esdata;
	InputGeneratorConfiguration generator;
	FETI4ILibraryConfiguration feti4ilibrary;

	HeatTransferConfiguration heat_transfer_2d;
	HeatTransferConfiguration heat_transfer_3d;
	StructuralMechanicsConfiguration structural_mechanics_2d;
	StructuralMechanicsConfiguration structural_mechanics_3d;

	OutputConfiguration output;

	DecomposerConfiguration decomposer;


	ECFConfiguration();
	ECFConfiguration(const std::string &file);
	ECFConfiguration(int *argc, char ***argv);
	bool fill(const std::string &file);
	bool fill(int *argc, char ***argv);

protected:
	void init();
};

}

#endif /* SRC_CONFIG_ECF_ECF_H_ */
