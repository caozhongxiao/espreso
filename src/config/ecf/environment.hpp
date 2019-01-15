
#ifndef SRC_CONFIG_ECF_ENVIRONMENT_HPP_
#define SRC_CONFIG_ECF_ENVIRONMENT_HPP_

#include "environment.h"
#include "config/description.h"

namespace espreso {

struct EnvironmentConfiguration: public Environment, public ECFDescription {

	EnvironmentConfiguration();

};

}



#endif /* SRC_CONFIG_ECF_ENVIRONMENT_HPP_ */
