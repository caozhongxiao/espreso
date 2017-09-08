
#include "output.h"
#include "../configuration.hpp"
#include "ecf.h"

espreso::ECFConfiguration* espreso::MonitorConfiguration::ecf = NULL;

espreso::MonitorConfiguration::MonitorConfiguration(const PHYSICS &physics)
: _physics(physics)
{
	REGISTER(region, ECFMetaData()
			.setdescription({ "Monitored region." })
			.setdatatype({ ECFDataType::REGION }));

	statistics = STATISTICS::AVG;
	REGISTER(statistics, ECFMetaData()
			.setdescription({ "Requested statistics." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("MIN").setdescription("Minimum."))
			.addoption(ECFOption().setname("MAX").setdescription("Maximum."))
			.addoption(ECFOption().setname("AVG").setdescription("Average."))
			.addoption(ECFOption().setname("NORM").setdescription("Norm (for testing purposes).").allowonly([] () { return false; })));

	auto temperature = [&] () {
		return _physics == PHYSICS::HEAT_TRANSFER_2D || _physics == PHYSICS::HEAT_TRANSFER_3D;
	};

	auto temperature3D = [&] () {
		return _physics == PHYSICS::HEAT_TRANSFER_3D;
	};

	auto structuralMechanics = [&] () {
		return _physics == PHYSICS::STRUCTURAL_MECHANICS_2D || _physics == PHYSICS::STRUCTURAL_MECHANICS_3D;
	};

	auto structuralMechanics3D = [&] () {
		return _physics == PHYSICS::STRUCTURAL_MECHANICS_3D;
	};

	REGISTER(property, ECFMetaData()
			.setdescription({ "Requested property." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("TEMPERATURE").setdescription("Minimum.").allowonly(temperature))

			.addoption(ECFOption().setname("FLUX").setdescription("Heat flux magnitude.").allowonly(temperature))
			.addoption(ECFOption().setname("FLUX_X").setdescription("Heat flux in x-direction.").allowonly(temperature))
			.addoption(ECFOption().setname("FLUX_Y").setdescription("Heat flux in y-direction.").allowonly(temperature))
			.addoption(ECFOption().setname("FLUX_Z").setdescription("Heat flux in z-direction.").allowonly(temperature3D))

			.addoption(ECFOption().setname("GRADIENT").setdescription("Heat gradient magnitude.").allowonly(temperature))
			.addoption(ECFOption().setname("GRADIENT_X").setdescription("Heat gradient in x-direction.").allowonly(temperature))
			.addoption(ECFOption().setname("GRADIENT_Y").setdescription("Heat gradient in y-direction.").allowonly(temperature))
			.addoption(ECFOption().setname("GRADIENT_Z").setdescription("Heat gradient in z-direction.").allowonly(temperature3D))

			.addoption(ECFOption().setname("DISPLACEMENT").setdescription("Displacement magnitude.").allowonly(structuralMechanics))
			.addoption(ECFOption().setname("DISPLACEMENT_X").setdescription("Displacement in x-direction.").allowonly(structuralMechanics))
			.addoption(ECFOption().setname("DISPLACEMENT_Y").setdescription("Displacement in y-direction.").allowonly(structuralMechanics))
			.addoption(ECFOption().setname("DISPLACEMENT_Z").setdescription("Displacement in z-direction.").allowonly(structuralMechanics3D))
	);
}

espreso::OutputConfiguration::OutputConfiguration(const PHYSICS &physics)
: _physics(physics)
{
	path = "results";
	REGISTER(path, ECFMetaData()
			.setdescription({ "Output path." })
			.setdatatype({ ECFDataType::STRING }));

	format = FORMAT::VTK_XML_ASCII;
	REGISTER(format, ECFMetaData()
			.setdescription({ "Output data format." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("VTK_LEGACY").setdescription("VTK legacy format."))
			.addoption(ECFOption().setname("VTK_XML_ASCII").setdescription("VTK XML ASCII format."))
			.addoption(ECFOption().setname("VTK_XML_BINARY").setdescription("VTK XML binary format."))
			.addoption(ECFOption().setname("ENSIGHT").setdescription("EnSight format.")));

	mode = MODE::THREAD;
	REGISTER(mode, ECFMetaData()
			.setdescription({ "ASYNC library mode." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("SYNC").setdescription("Output is synchronized."))
			.addoption(ECFOption().setname("THREAD").setdescription("Output is done by asynchronous thread."))
			.addoption(ECFOption().setname("MPI").setdescription("Output is done by separated MPI processes.")));

	output_node_group_size = 7;
	REGISTER(output_node_group_size, ECFMetaData()
			.setdescription({ "Number of MPI processes that send output data to the same storing node." })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	solution = true;
	subsolution = settings = FETI_data = catalyst = false;
	catalyst_sleep_time = 0;
	REGISTER(solution, ECFMetaData()
			.setdescription({ "Store solution." })
			.setdatatype({ ECFDataType::BOOL }));
	REGISTER(subsolution, ECFMetaData()
			.setdescription({ "Store subsolution of non-linear solvers." })
			.setdatatype({ ECFDataType::BOOL }));
	REGISTER(settings, ECFMetaData()
			.setdescription({ "Store settings." })
			.setdatatype({ ECFDataType::BOOL }));
	REGISTER(FETI_data, ECFMetaData()
			.setdescription({ "Store FETI related data." })
			.setdatatype({ ECFDataType::BOOL }));
	REGISTER(catalyst, ECFMetaData()
			.setdescription({ "In-situ visualization by Catalyst." })
			.setdatatype({ ECFDataType::BOOL }));
	REGISTER(catalyst_sleep_time, ECFMetaData()
			.setdescription({ "The sleep time between each time steps when catalyst is used." })
			.setdatatype({ ECFDataType::NONNEGATIVE_INTEGER }));

	collected = separate_bodies = separate_materials = false;
	REGISTER(collected, ECFMetaData()
			.setdescription({ "Collect the results to one file." })
			.setdatatype({ ECFDataType::BOOL }));
	REGISTER(separate_bodies, ECFMetaData()
			.setdescription({ "Separate bodies to regions." })
			.setdatatype({ ECFDataType::BOOL }));
	REGISTER(separate_materials, ECFMetaData()
			.setdescription({ "Separate materials to regions." })
			.setdatatype({ ECFDataType::BOOL }));

	domain_shrink_ratio = .95;
	cluster_shrink_ratio = .9;
	REGISTER(domain_shrink_ratio, ECFMetaData()
			.setdescription({ "Domain shrink ratio." })
			.setdatatype({ ECFDataType::FLOAT }));
	REGISTER(cluster_shrink_ratio, ECFMetaData()
			.setdescription({ "Cluster shrink ratio." })
			.setdatatype({ ECFDataType::FLOAT }));

	REGISTER(
			monitoring,
			ECFMetaData()
				.setdescription({ "List of monitored data.", "Column index of *.emr file." })
				.setdatatype({ ECFDataType::POSITIVE_INTEGER })
				.setpattern({ "1" }),
			_physics);
}



