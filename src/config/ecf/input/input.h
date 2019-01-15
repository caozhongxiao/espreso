
#ifndef SRC_CONFIG_ECF_INPUT_INPUT_H_
#define SRC_CONFIG_ECF_INPUT_INPUT_H_

#include "config/description.h"
#include "config/ecf/processesreduction.h"

#include <string>

namespace espreso {

enum class INPUT_FORMAT {
	WORKBENCH,
	OPENFOAM,
	ESDATA,
	GENERATOR
};


struct InputConfiguration: public ProcessesReduction, public ECFDescription {

	std::string path;

	bool keep_material_sets;
	bool convert_database;

	double scale_factor;

	InputConfiguration();
};

}



#endif /* SRC_CONFIG_ECF_INPUT_INPUT_H_ */
