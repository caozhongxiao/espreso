
#ifndef SRC_CONFIG_ECF_ROOT_H_
#define SRC_CONFIG_ECF_ROOT_H_

#include "environment.h"
#include "output.h"
#include "decomposer.h"
#include "decomposition.h"

#include "pythontestgenerator.h"

#include "input/input.h"
#include "input/generator.h"
#include "input/feti4ilibrary.h"

#include "meshmorphing.h"

#include "physics/physics.h"
#include "physics/structuralmechanics.h"

#include "../reader/reader.h"
#include "physics/heattransfer.h"

namespace espreso {

struct ECFRoot: public ECFObject {

	ECFObject* getInput() { return const_cast<ECFObject*>(_getInput()); }
	const ECFObject* getInput() const { return _getInput(); }

	PhysicsConfiguration* getPhysics() { return const_cast<PhysicsConfiguration*>(_getPhysics()); }
	const PhysicsConfiguration* getPhysics() const { return _getPhysics(); }

	// Environment has to be created first!
	Environment environment;

	PythonTestGenerator python_test_generator;

	std::map<size_t, std::string> default_args;
	std::map<std::string, std::string> variables;

	INPUT_FORMAT input;
	PHYSICS physics;

	DecompositionConfiguration decomposition;

	InputConfiguration workbench, openfoam, esdata;
	InputGeneratorConfiguration generator;
	FETI4ILibraryConfiguration feti4ilibrary;

	MeshMorphing mesh_morphing;

	HeatTransferConfiguration heat_transfer_2d;
	HeatTransferConfiguration heat_transfer_3d;
	StructuralMechanicsConfiguration structural_mechanics_2d;
	StructuralMechanicsConfiguration structural_mechanics_3d;

	OutputConfiguration output;

	DecomposerConfiguration decomposer;

	ECFRoot();
	ECFRoot(const std::string &file);
	ECFRoot(int *argc, char ***argv);
	bool fill(const std::string &file);
	bool fill(int *argc, char ***argv);

protected:
	void init();

	const ECFObject* _getInput() const;
	const PhysicsConfiguration* _getPhysics() const;
};

}

#endif /* SRC_CONFIG_ECF_ROOT_H_ */
