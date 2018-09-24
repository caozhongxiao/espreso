
#ifndef SRC_CONFIG_ECF_INPUT_INPUT_H_
#define SRC_CONFIG_ECF_INPUT_INPUT_H_

#include "loaderconfiguration.h"
#include "../../configuration.h"

namespace espreso {

enum class INPUT_FORMAT {
	WORKBENCH,
	OPENFOAM,
	ESDATA,
	GENERATOR
};


struct InputConfiguration: public LoaderConfiguration, public ECFObject {

	std::string path;

	bool keep_material_sets;
	bool convert_database;

	double scale_factor;

	InputConfiguration();
};

}



#endif /* SRC_CONFIG_ECF_INPUT_INPUT_H_ */
