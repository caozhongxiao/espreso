
#ifndef SRC_CONFIG_ECF_DECOMPOSITION_H_
#define SRC_CONFIG_ECF_DECOMPOSITION_H_

#include "../configuration.h"

namespace espreso {

struct DecompositionConfiguration: public ECFObject {

	std::string path;

	int mpi_procs, domains;

	bool balance_clusters;

	bool separate_materials, separate_regions, separate_etypes;

	DecompositionConfiguration();
};

}



#endif /* SRC_CONFIG_ECF_DECOMPOSITION_H_ */
