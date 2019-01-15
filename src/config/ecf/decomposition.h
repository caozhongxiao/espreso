
#ifndef SRC_CONFIG_ECF_DECOMPOSITION_H_
#define SRC_CONFIG_ECF_DECOMPOSITION_H_

#include "config/description.h"
#include "processesreduction.h"

#include <string>

namespace espreso {

struct METISConfiguration: public ProcessesReduction, public ECFDescription {

	enum class OBJECTIVE_TYPE {
		VOLUME, EDGECUT
	};

	OBJECTIVE_TYPE objective_type;
	bool refinement;

	METISConfiguration();
};

struct DecompositionConfiguration: public ECFDescription {

	std::string path;

	int mpi_procs, domains;

	bool balance_clusters;

	bool separate_materials, separate_regions, separate_etypes;

	METISConfiguration metis_options;

	DecompositionConfiguration();
};

}



#endif /* SRC_CONFIG_ECF_DECOMPOSITION_H_ */
