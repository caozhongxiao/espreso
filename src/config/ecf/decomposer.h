
#ifndef SRC_CONFIG_ECF_DECOMPOSER_H_
#define SRC_CONFIG_ECF_DECOMPOSER_H_

#include "../configuration.h"

namespace espreso {

struct DecomposerConfiguration: public ECFObject {

	std::string prefix, parts;

	DecomposerConfiguration();
};

}


#endif /* SRC_CONFIG_ECF_DECOMPOSER_H_ */
