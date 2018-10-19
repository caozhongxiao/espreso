
#ifndef SRC_CONFIG_ECF_OUTPUT_H_
#define SRC_CONFIG_ECF_OUTPUT_H_

#include "../configuration.h"
#include "physics/physics.h"

#include <iostream>

namespace espreso {

struct ECFRoot;

struct MonitorConfiguration: public ECFObject {

	enum class STATISTICS {
		MIN,
		MAX,
		AVG,
		NORM
	};

	std::string region;
	STATISTICS statistics;
	std::string property;

	MonitorConfiguration(const PHYSICS &physics);
protected:
	const PHYSICS &_physics;
};

struct ResultsSelectionConfiguration: public ECFObject {

	bool temperature, gradient, flux, phase, latent_heat;
	bool displacement;

	void basic();
	void all();

	ResultsSelectionConfiguration(const PHYSICS &physics);
protected:
	const PHYSICS &_physics;
};

struct OutputConfiguration: public ECFObject {

	enum class FORMAT {
		VTK_LEGACY = 0,
		VTK_XML_ASCII = 1,
		VTK_XML_BINARY = 2,
		ENSIGHT = 3,
		STL_SURFACE = 4
	};

	enum class MODE {
		SYNC,
		THREAD,
//		MPI,
	};

	enum class STORE_FREQUENCY {
		NEVER,
		EVERY_TIMESTEP,
		EVERY_NTH_TIMESTEP,
		LAST_TIMESTEP,
		DEBUG
	};

	enum class STORE_RESULTS {
		BASIC,
		ALL,
		USER
	};

	FORMAT format;
	MODE mode;

//	size_t output_node_group_size;

	std::string path;

	STORE_FREQUENCY results_store_frequency, monitors_store_frequency;
	size_t results_nth_stepping, monitors_nth_stepping;

	STORE_RESULTS store_results;
	ResultsSelectionConfiguration results_selection;

	bool settings, debug, catalyst;
	size_t catalyst_sleep_time;

	bool collected;

	std::map<size_t, MonitorConfiguration> monitoring;

	OutputConfiguration(const PHYSICS &physics);

protected:
	const PHYSICS &_physics;
};

}



#endif /* SRC_CONFIG_ECF_OUTPUT_H_ */
