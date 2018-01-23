
#ifndef SRC_CONFIG_ECF_INPUT_INPUT_H_
#define SRC_CONFIG_ECF_INPUT_INPUT_H_

#include "../../configuration.h"

namespace espreso {

enum class INPUT_FORMAT {
	WORKBENCH,
	OPENFOAM,
	ESDATA,
	GENERATOR
};

struct InputConfiguration: public ECFObject {

	std::string path;

	bool compress_numbers;

	InputConfiguration();
};

}



#endif /* SRC_CONFIG_ECF_INPUT_INPUT_H_ */
