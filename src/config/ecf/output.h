
#ifndef SRC_CONFIG_ECF_OUTPUT_H_
#define SRC_CONFIG_ECF_OUTPUT_H_

#include "config/description.h"
#include "physics/physics.h"

#include "physics/heattransfer.h"
#include "physics/structuralmechanics.h"

namespace espreso {

struct ECFRoot;

struct MonitorConfiguration: public ECFDescription {

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

struct ResultsSelectionConfiguration: public HeatTransferOutputSettings, public StructuralMechanicsOutputSettings, public ECFDescription {

	void basic()
	{
		HeatTransferOutputSettings::basic();
		StructuralMechanicsOutputSettings::basic();
	}

	void all()
	{
		HeatTransferOutputSettings::all();
		StructuralMechanicsOutputSettings::all();
	}

	ResultsSelectionConfiguration(const PHYSICS &physics);
protected:
	const PHYSICS &_physics;
};

struct OutputConfiguration: public ECFDescription {

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
		LAST_TIMESTEP
	};

	enum class STORE_RESULTS {
		BASIC,
		ALL,
		USER
	};

	std::string log_dir;

	size_t verbose_level;
	size_t measure_level;

	bool print_matrices;

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
