
#ifndef SRC_CONFIG_ECF_OUTPUT_H_
#define SRC_CONFIG_ECF_OUTPUT_H_

#include "../configuration.h"

namespace espreso {

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

	// SUBMULTIMAP(std::string, std::string, monitoring, "Results statistics in some regions. OPERATION = { AVERAGE, MIN, MAX }", "REGION", "<OPERATION> <VARIABLE>");

	OutputConfiguration();
};

}



#endif /* SRC_CONFIG_ECF_OUTPUT_H_ */
