
#include "output.h"
#include "../configuration.hpp"
#include "root.h"

espreso::MonitorConfiguration::MonitorConfiguration(const PHYSICS &physics)
: _physics(physics)
{
	REGISTER(region, ECFMetaData()
            .setdescription({ "Region" })
			.setdatatype({ ECFDataType::REGION }));

	statistics = STATISTICS::AVG;
	REGISTER(statistics, ECFMetaData()
            .setdescription({ "Statistics" })
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
            .setdescription({ "Result" })
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

espreso::ResultsSelectionConfiguration::ResultsSelectionConfiguration(const PHYSICS &physics)
: _physics(physics)
{
	basic();

	auto thermal = [&] () {
		return _physics == PHYSICS::HEAT_TRANSFER_2D || _physics == PHYSICS::HEAT_TRANSFER_3D;
	};

	auto structuralMechanics = [&] () {
		return _physics == PHYSICS::STRUCTURAL_MECHANICS_2D || _physics == PHYSICS::STRUCTURAL_MECHANICS_3D;
	};

	REGISTER(temperature, ECFMetaData()
			.setdescription({ "Temperature." })
			.setdatatype({ ECFDataType::BOOL })
			.allowonly(thermal));
	REGISTER(gradient, ECFMetaData()
			.setdescription({ "Temperature gradient." })
			.setdatatype({ ECFDataType::BOOL })
			.allowonly(thermal));
	REGISTER(flux, ECFMetaData()
			.setdescription({ "Temperature flux." })
			.setdatatype({ ECFDataType::BOOL })
			.allowonly(thermal));

	REGISTER(phase, ECFMetaData()
			.setdescription({ "Material phase." })
			.setdatatype({ ECFDataType::BOOL })
			.allowonly(thermal));
	REGISTER(latent_heat, ECFMetaData()
			.setdescription({ "Latent heat." })
			.setdatatype({ ECFDataType::BOOL })
			.allowonly(thermal));

	REGISTER(displacement, ECFMetaData()
			.setdescription({ "Displacement." })
			.setdatatype({ ECFDataType::BOOL })
			.allowonly(structuralMechanics));
}

void espreso::ResultsSelectionConfiguration::basic()
{
	temperature = true;
	gradient = flux = phase = latent_heat = false;
	displacement = true;
}

void espreso::ResultsSelectionConfiguration::all()
{
	temperature = gradient = flux = phase = latent_heat = true;
	displacement = true;
}

espreso::OutputConfiguration::OutputConfiguration(const PHYSICS &physics)
: results_selection(physics), _physics(physics)
{
	path = "results";
	REGISTER(path, ECFMetaData()
            .setdescription({ "Path" })
			.setdatatype({ ECFDataType::STRING }));

	format = FORMAT::VTK_XML_ASCII;
	REGISTER(format, ECFMetaData()
            .setdescription({ "Format" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("VTK_LEGACY").setdescription("VTK legacy format."))
			.addoption(ECFOption().setname("VTK_XML_ASCII").setdescription("VTK XML ASCII format."))
			.addoption(ECFOption().setname("VTK_XML_BINARY").setdescription("VTK XML binary format."))
			.addoption(ECFOption().setname("ENSIGHT").setdescription("EnSight format.")));

	mode = MODE::THREAD;
	REGISTER(mode, ECFMetaData()
            .setdescription({ "ASYNC library mode" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("SYNC").setdescription("Output is synchronized."))
			.addoption(ECFOption().setname("THREAD").setdescription("Output is done by asynchronous thread.")));
			//.addoption(ECFOption().setname("MPI").setdescription("Output is done by separated MPI processes.")));

//	output_node_group_size = 7;
//	REGISTER(output_node_group_size, ECFMetaData()
//			.setdescription({ "Number of MPI processes that send output data to the same storing node." })
//			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	addSeparator();

	results_store_frequency = monitors_store_frequency = STORE_FREQUENCY::EVERY_TIMESTEP;
	results_nth_stepping = monitors_nth_stepping = 10;
	REGISTER(results_store_frequency, ECFMetaData()
            .setdescription({ "Results store frequency" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("NEVER").setdescription("Storing is turned off."))
			.addoption(ECFOption().setname("EVERY_TIMESTEP").setdescription("Results are stored after each time step."))
			.addoption(ECFOption().setname("EVERY_NTH_TIMESTEP").setdescription("Results are stored after each nth time step."))
			.addoption(ECFOption().setname("LAST_TIMESTEP").setdescription("Only last results are stored."))
			.addoption(ECFOption().setname("DEBUG").setdescription("Storing also iteration results.")));
	REGISTER(monitors_store_frequency, ECFMetaData()
            .setdescription({ "Monitoring store frequency" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("NEVER").setdescription("Storing is turned off."))
			.addoption(ECFOption().setname("EVERY_TIMESTEP").setdescription("Monitors are stored after each time step."))
			.addoption(ECFOption().setname("EVERY_NTH_TIMESTEP").setdescription("Monitors are stored after each nth time step."))
			.addoption(ECFOption().setname("LAST_TIMESTEP").setdescription("Only last monitors are stored.")));

	REGISTER(results_nth_stepping, ECFMetaData()
            .setdescription({ "Write results" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER })
			.allowonly([&] () { return results_store_frequency == STORE_FREQUENCY::EVERY_NTH_TIMESTEP; }));
	REGISTER(monitors_nth_stepping, ECFMetaData()
            .setdescription({ "Monitors store stepping" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER})
			.allowonly([&] () { return monitors_store_frequency == STORE_FREQUENCY::EVERY_NTH_TIMESTEP; }));

	addSpace();

	store_results = STORE_RESULTS::BASIC;
	REGISTER(store_results, ECFMetaData()
            .setdescription({ "Stored properties" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("BASIC").setdescription("Basic properties."))
			.addoption(ECFOption().setname("ALL").setdescription("All properties."))
			.addoption(ECFOption().setname("USER").setdescription("User defined properties.")))
	->addListener(ECFParameter::Event::VALUE_SET, [&] (const std::string &value) {
		switch (store_results) {
		case STORE_RESULTS::BASIC:
			results_selection.basic();
			break;
		case STORE_RESULTS::ALL:
			results_selection.all();
			break;
		case STORE_RESULTS::USER:
			break;
		}
	});

	REGISTER(results_selection, ECFMetaData()
            .setdescription({ "Properties selection" })
			.allowonly([&] () { return store_results == STORE_RESULTS::USER; }));

	addSeparator();

	settings = debug = false;
	catalyst = false;
	catalyst_sleep_time = 0;
	REGISTER(settings, ECFMetaData()
            .setdescription({ "Store settings" })
			.setdatatype({ ECFDataType::BOOL }));
	REGISTER(debug, ECFMetaData()
            .setdescription({ "Store FETI related data" })
			.setdatatype({ ECFDataType::BOOL }));
	REGISTER(catalyst, ECFMetaData()
            .setdescription({ "In-situ visualization by Catalyst" })
			.setdatatype({ ECFDataType::BOOL }));
	REGISTER(catalyst_sleep_time, ECFMetaData()
            .setdescription({ "The sleep time between each time steps when catalyst is used" })
			.setdatatype({ ECFDataType::NONNEGATIVE_INTEGER }));

	addSpace();

	collected = true;
	REGISTER(collected, ECFMetaData()
            .setdescription({ "Collected results" })
			.setdatatype({ ECFDataType::BOOL }));

	addSpace();

	REGISTER(
			monitoring,
			ECFMetaData()
                .setdescription({ "List of monitors", "Monitor" })
				.setdatatype({ ECFDataType::POSITIVE_INTEGER })
				.setpattern({ "1" }),
			_physics);
}



