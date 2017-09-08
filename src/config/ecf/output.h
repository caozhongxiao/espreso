
#ifndef SRC_CONFIG_ECF_OUTPUT_H_
#define SRC_CONFIG_ECF_OUTPUT_H_

#include "../configuration.h"
#include "physics/physics.h"

#include <iostream>

namespace espreso {

struct ECFConfiguration;

struct MonitorConfiguration: public ECFObject {

	static ECFConfiguration *ecf;

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

struct OutputConfiguration: public ECFObject {

	enum class FORMAT {
		VTK_LEGACY = 0,
		VTK_XML_ASCII = 1,
		VTK_XML_BINARY = 2,
		ENSIGHT = 3
	};

	enum class MODE {
		SYNC,
		THREAD,
		MPI,
	};

	FORMAT format;
	MODE mode;

	size_t output_node_group_size;

	std::string path;

	bool solution, subsolution, settings, FETI_data, catalyst;
	size_t catalyst_sleep_time;

	bool collected, separate_bodies, separate_materials;
	double domain_shrink_ratio, cluster_shrink_ratio;

	std::map<size_t, MonitorConfiguration> monitoring;

	OutputConfiguration(const PHYSICS &physics);

protected:
	const PHYSICS &_physics;
};

}



#endif /* SRC_CONFIG_ECF_OUTPUT_H_ */
