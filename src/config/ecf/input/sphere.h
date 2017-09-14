
#ifndef SRC_CONFIG_ECF_INPUT_SPHERE_H_
#define SRC_CONFIG_ECF_INPUT_SPHERE_H_

#include "generatorelements.h"
#include "../../configuration.h"

namespace espreso {

struct SphereGeneratorConfiguration: public ECFObject {

	GENERATOR_ELEMENT_TYPE element_type;

	double inner_radius, outer_radius;
	size_t clusters, layers;
	size_t domains_x, domains_y, domains_z;
	size_t elements_x, elements_y, elements_z;

	bool uniform_decomposition;

	std::map<std::string, std::string> nodes, edges, faces, elements;

	SphereGeneratorConfiguration();
};

}



#endif /* SRC_CONFIG_ECF_INPUT_SPHERE_H_ */
