
#ifndef SRC_CONFIG_ECF_INPUT_INPUT_H_
#define SRC_CONFIG_ECF_INPUT_INPUT_H_

#include "../../configuration.h"
#include "../processesreduction.h"

namespace espreso {

enum class INPUT_FORMAT {
	WORKBENCH,
	OPENFOAM,
	ABAQUS,
	ESDATA,
	GENERATOR
};


struct InputConfiguration: public ProcessesReduction, public ECFObject {

	std::string path;

	bool keep_material_sets;
	bool convert_database;

	double scale_factor;

	InputConfiguration();
};

}



#endif /* SRC_CONFIG_ECF_INPUT_INPUT_H_ */
