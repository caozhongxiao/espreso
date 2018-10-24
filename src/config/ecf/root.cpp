
#include "root.h"

#include "../configuration.hpp"

#include "../../basis/logging/logging.h"

using namespace espreso;

const ECFObject* ECFRoot::_getInput() const
{
	switch (input) {
	case INPUT_FORMAT::GENERATOR:
		return &generator;
	case INPUT_FORMAT::WORKBENCH:
		return &workbench;
	case INPUT_FORMAT::OPENFOAM:
		return &openfoam;
	case INPUT_FORMAT::ESDATA:
		return &esdata;
	default:
		ESINFO(GLOBAL_ERROR) << "Request for unknown input.";
		return NULL;
	}
}

const PhysicsConfiguration* ECFRoot::_getPhysics() const
{
	switch (physics) {
	case PHYSICS::HEAT_TRANSFER_2D:
		return &heat_transfer_2d;
	case PHYSICS::HEAT_TRANSFER_3D:
		return &heat_transfer_3d;
	case PHYSICS::STRUCTURAL_MECHANICS_2D:
		return &structural_mechanics_2d;
	case PHYSICS::STRUCTURAL_MECHANICS_3D:
		return &structural_mechanics_3d;
	default:
		ESINFO(GLOBAL_ERROR) << "Request for unknown physics.";
		return NULL;
	}
}


void ECFRoot::init()
{
	name = "root";

	REGISTER(python_test_generator, ECFMetaData()
			.setdescription({ "Description of Python test generator (run python tests/generate.py PATH)." }));

	REGISTER(default_args, ECFMetaData()
		.setdescription({ "The index of the argument.", "The argument value." })
		.setdatatype({ ECFDataType::NONNEGATIVE_INTEGER, ECFDataType::STRING })
		.setpattern({ "0", "VALUE" }));

	REGISTER(variables, ECFMetaData()
			.setdescription({ "A name of variable usable in *.ecf file.", "A value of the variable." })
			.setdatatype({ ECFDataType::STRING, ECFDataType::STRING })
			.setpattern({ "MY_VARIABLE", "VALUE" }));

	addSpace();

	REGISTER(decomposition, ECFMetaData()
			.setdescription({ "Domains decomposition settings." }));

	input = INPUT_FORMAT::GENERATOR;

	REGISTER(input, ECFMetaData()
			.setdescription({ "An input data type." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("WORKBENCH").setdescription("Ansys WorkBench format."))
			.addoption(ECFOption().setname("OPENFOAM").setdescription("OpenFOAM format."))
			.addoption(ECFOption().setname("ESDATA").setdescription("ESPRESO internal binary format."))
			.addoption(ECFOption().setname("GENERATOR").setdescription("ESPRESO internal generator.")));

	physics = PHYSICS::HEAT_TRANSFER_3D;
	REGISTER(physics, ECFMetaData()
			.setdescription({ "Physics" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("HEAT_TRANSFER_2D").setdescription("Heat transfer 2D."))
			.addoption(ECFOption().setname("HEAT_TRANSFER_3D").setdescription("Heat transfer 3D."))
			.addoption(ECFOption().setname("STRUCTURAL_MECHANICS_2D").setdescription("Structural mechanics 2D."))
			.addoption(ECFOption().setname("STRUCTURAL_MECHANICS_3D").setdescription("Structural mechanics 3D.")));

	REGISTER(workbench, ECFMetaData()
			.setdescription({ "Description of Ansys WorkBench format." })
			.allowonly([&] () { return input == INPUT_FORMAT::WORKBENCH; }));
	REGISTER(openfoam, ECFMetaData()
			.setdescription({ "Description of OpenFoam format." })
			.allowonly([&] () { return input == INPUT_FORMAT::OPENFOAM; }));
	REGISTER(esdata, ECFMetaData()
			.setdescription({ "Description of ESDATA format." })
			.allowonly([&] () { return input == INPUT_FORMAT::ESDATA; }));
	REGISTER(generator, ECFMetaData()
			.setdescription({ "Description of ESPRESO generator." })
			.allowonly([&] () { return input == INPUT_FORMAT::GENERATOR; }));

	REGISTER(feti4ilibrary, ECFMetaData()
			.setdescription({ "Settings for FETI4I library." }));

	REGISTER(mesh_morphing, ECFMetaData()
			.setdescription({ "Settings for mesh morphing." }));

	REGISTER(heat_transfer_2d, ECFMetaData()
			.setdescription({ "Heat transfer 2D" })
			.allowonly([&] () { return physics == PHYSICS::HEAT_TRANSFER_2D; }));
	REGISTER(heat_transfer_3d, ECFMetaData()
			.setdescription({ "Heat transfer 3D" })
			.allowonly([&] () { return physics == PHYSICS::HEAT_TRANSFER_3D; }));
	REGISTER(structural_mechanics_2d, ECFMetaData()
			.setdescription({ "Structural mechanics 2D" })
			.allowonly([&] () { return physics == PHYSICS::STRUCTURAL_MECHANICS_2D; }));
	REGISTER(structural_mechanics_3d, ECFMetaData()
			.setdescription({ "Structural mechanics 3D" })
			.allowonly([&] () { return physics == PHYSICS::STRUCTURAL_MECHANICS_3D; }));

	registerParameter("env", environment, ECFMetaData()
			.setdescription({ "Environment related settings." }));

	REGISTER(output, ECFMetaData()
			.setdescription({ "Output configurations." }));
}

ECFRoot::ECFRoot()
: mesh_morphing(this),
  heat_transfer_2d(DIMENSION::D2),
  heat_transfer_3d(DIMENSION::D3),
  structural_mechanics_2d(DIMENSION::D2),
  structural_mechanics_3d(DIMENSION::D3),
  output(physics)
{
	init();
}

ECFRoot::ECFRoot(const std::string &file)
: ECFRoot()
{
	if (!fill(file)) {
		ESINFO(GLOBAL_ERROR) << "Cannot read ECF file '" << file << "'.";
	}
}

ECFRoot::ECFRoot(int *argc, char ***argv)
: ECFRoot()
{
	if (!fill(argc, argv)) {
		ESINFO(GLOBAL_ERROR)
				<< "Cannot read ECF file '" << ECFReader::configurationFile << "'.\n"
				<< "Use default 'espreso.ecf' or set an arbitrary by 'espreso -c $ecfpath'.";
	}
}

bool ECFRoot::fill(const std::string &file)
{
	if (ECFReader::read(*this, file, this->default_args, this->variables).hadValidECF) {
		ECFReader::set(this->environment, this->output);
		return true;
	}
	return false;
}

bool ECFRoot::fill(int *argc, char ***argv)
{
	if (ECFReader::read(*this, argc, argv, this->default_args, this->variables).hadValidECF) {
		ECFReader::set(this->environment, this->output);
		return true;
	}
	return false;
}



