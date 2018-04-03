
#ifndef SRC_OUTPUT_RESULT_VISUALIZATION_VISUALIZATION_H_
#define SRC_OUTPUT_RESULT_VISUALIZATION_VISUALIZATION_H_

#include "../resultstore.h"

namespace espreso {

class Visualization: public ResultStoreBase {

public:
	static bool storeStep(const OutputConfiguration &configuration, const Step &step);

	Visualization(const Mesh &mesh, const OutputConfiguration &configuration);

protected:
	const OutputConfiguration &_configuration;
};

}



#endif /* SRC_OUTPUT_RESULT_VISUALIZATION_VISUALIZATION_H_ */
