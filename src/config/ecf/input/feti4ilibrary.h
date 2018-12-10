
#ifndef SRC_CONFIG_ECF_INPUT_FETI4ILIBRARY_H_
#define SRC_CONFIG_ECF_INPUT_FETI4ILIBRARY_H_

#include "../linearsolver/feti.h"

namespace espreso {

struct FETI4ILibraryConfiguration: public ECFObject {

	size_t domains;
	FETISolverConfiguration solver;

	FETI4ILibraryConfiguration();
};

}



#endif /* SRC_CONFIG_ECF_INPUT_FETI4ILIBRARY_H_ */
