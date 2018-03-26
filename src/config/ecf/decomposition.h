
#ifndef SRC_CONFIG_ECF_DECOMPOSITION_H_
#define SRC_CONFIG_ECF_DECOMPOSITION_H_

#include "../configuration.h"

namespace espreso {

struct METISConfiguration: public ECFObject {

	enum class OBJECTIVE_TYPE {
		VOLUME, EDGECUT
	};

	OBJECTIVE_TYPE objective_type;
	bool adaptive_refinement;

	METISConfiguration();
};

struct DecompositionConfiguration: public ECFObject {

	std::string path;

	int mpi_procs, domains;

	bool balance_clusters;

	bool separate_materials, separate_regions, separate_etypes;

	METISConfiguration metis_options;

	DecompositionConfiguration();
};

}



#endif /* SRC_CONFIG_ECF_DECOMPOSITION_H_ */
