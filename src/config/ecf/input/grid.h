
#ifndef SRC_CONFIG_ECF_INPUT_GRID_H_
#define SRC_CONFIG_ECF_INPUT_GRID_H_

#include "generatorelements.h"
#include "../../configuration.h"

namespace espreso {

struct GridGeneratorConfiguration: public ECFObject {

	GENERATOR_ELEMENT_TYPE element_type;

	double start_x, start_y, start_z;
	double length_x, length_y, length_z;

	std::string projection_x, projection_y, projection_z;
	std::string rotation_x, rotation_y, rotation_z;

	size_t blocks_x, blocks_y, blocks_z;
	size_t clusters_x, clusters_y, clusters_z;
	size_t domains_x, domains_y, domains_z;
	size_t elements_x, elements_y, elements_z;

	bool random_partition, uniform_decomposition;

	std::map<size_t, bool> blocks;
	std::map<size_t, size_t> noncontinuous;

	std::map<std::string, std::string> nodes, edges, faces, elements;

	size_t chessboard_size;
	size_t nonuniform_nparts;

	GridGeneratorConfiguration();
};

}



#endif /* SRC_CONFIG_ECF_INPUT_GRID_H_ */
