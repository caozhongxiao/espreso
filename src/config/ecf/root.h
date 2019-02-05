
#ifndef SRC_CONFIG_ECF_ROOT_H_
#define SRC_CONFIG_ECF_ROOT_H_

#include "output.h"
#include "decomposition.h"

#include "pythontestgenerator.h"

#include "input/input.h"
#include "input/generator.h"
#include "input/feti4ilibrary.h"

#include "meshmorphing.h"

#include "physics/physics.h"
#include "physics/heattransfer.h"
#include "physics/structuralmechanics.h"

#include "config/reader/reader.h"


namespace espreso {

struct ECFRoot: public ECFDescription {

	ECFDescription* getInput() { return const_cast<ECFDescription*>(_getInput()); }
	const ECFDescription* getInput() const { return _getInput(); }

	PhysicsConfiguration* getPhysics() { return const_cast<PhysicsConfiguration*>(_getPhysics()); }
	const PhysicsConfiguration* getPhysics() const { return _getPhysics(); }

	PythonTestGenerator python_test_generator;

	std::map<size_t, std::string> default_args;
	std::map<std::string, std::string> variables;

	INPUT_FORMAT input;
	PHYSICS physics;

	DecompositionConfiguration decomposition;

	InputConfiguration workbench, openfoam, abaqus, esdata;
	InputGeneratorConfiguration generator;
	FETI4ILibraryConfiguration feti4ilibrary;

	MeshMorphing mesh_morphing;

	HeatTransferConfiguration heat_transfer_2d;
	HeatTransferConfiguration heat_transfer_3d;
	StructuralMechanicsConfiguration structural_mechanics_2d;
	StructuralMechanicsConfiguration structural_mechanics_3d;

	OutputConfiguration output;

	ECFRoot();
	ECFRoot(const std::string &file);
	ECFRoot(int *argc, char ***argv);
	bool fill(const std::string &file);
	bool fill(int *argc, char ***argv);

protected:
	void init();

	const ECFDescription* _getInput() const;
	const PhysicsConfiguration* _getPhysics() const;
};

}

#endif /* SRC_CONFIG_ECF_ROOT_H_ */
