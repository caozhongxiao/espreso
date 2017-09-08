
#include "ecf.h"

#include "../configuration.hpp"

void espreso::ECFConfiguration::init()
{
	name = "root";

	REGISTER(default_args, ECFMetaData()
		.setdescription({ "The index of the argument.", "The argument value." })
		.setdatatype({ ECFDataType::NONNEGATIVE_INTEGER, ECFDataType::STRING })
		.setpattern({ "0", "VALUE" }));

	REGISTER(variables, ECFMetaData()
			.setdescription({ "A name of variable usable in *.ecf file.", "A value of the variable." })
			.setdatatype({ ECFDataType::STRING, ECFDataType::STRING })
			.setpattern({ "MY_VARIABLE", "VALUE" }));

	addSpace();

	input = INPUT_FORMAT::GENERATOR;

	REGISTER(input, ECFMetaData()
			.setdescription({ "An input data type." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("WORKBENCH").setdescription("Ansys WorkBench format."))
			.addoption(ECFOption().setname("OPENFOAM").setdescription("OpenFOAM format."))
			.addoption(ECFOption().setname("ESDATA").setdescription("ESPRESO internal binary format."))
			.addoption(ECFOption().setname("GENERATOR").setdescription("ESPRESO internal generator.")));

	physics = PHYSICS::ADVECTION_DIFFUSION_2D;
	REGISTER(physics, ECFMetaData()
			.setdescription({ "A selected physics." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("ADVECTION_DIFFUSION_2D").setdescription("Advection diffusion 2D."))
			.addoption(ECFOption().setname("ADVECTION_DIFFUSION_3D").setdescription("Advection diffusion 3D."))
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

	REGISTER(advection_diffusion_2d, ECFMetaData()
			.setdescription({ "Advection diffusion 2D settings." })
			.allowonly([&] () { return physics == PHYSICS::ADVECTION_DIFFUSION_2D; }));
	REGISTER(advection_diffusion_3d, ECFMetaData()
			.setdescription({ "Advection diffusion 3D settings." })
			.allowonly([&] () { return physics == PHYSICS::ADVECTION_DIFFUSION_3D; }));
	REGISTER(structural_mechanics_2d, ECFMetaData()
			.setdescription({ "Structural mechanics 2D settings." })
			.allowonly([&] () { return physics == PHYSICS::STRUCTURAL_MECHANICS_2D; }));
	REGISTER(structural_mechanics_3d, ECFMetaData()
			.setdescription({ "Structural mechanics 3D settings." })
			.allowonly([&] () { return physics == PHYSICS::STRUCTURAL_MECHANICS_3D; }));

	registerParameter("env", environment, ECFMetaData()
			.setdescription({ "Environment related settings." }));

	REGISTER(output, ECFMetaData()
			.setdescription({ "Output configurations." }));


	REGISTER(decomposer, ECFMetaData()
			.setdescription({ "Configuration of ESPRESO decomposer." }));
}

espreso::ECFConfiguration::ECFConfiguration()
: output(physics)
{
	init();
}

espreso::ECFConfiguration::ECFConfiguration(const std::string &file)
: output(physics)
{
	init();
	ECFReader::read(*this, file, this->default_args, this->variables);
	ECFReader::set(this->environment, this->output);
}

espreso::ECFConfiguration::ECFConfiguration(int *argc, char ***argv)
: output(physics)
{
	init();
	ECFReader::read(*this, argc, argv, this->default_args, this->variables);
	ECFReader::set(this->environment, this->output);
}



