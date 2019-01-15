
#ifndef SRC_OUTPUT_RESULT_VISUALIZATION_VISUALIZATION_H_
#define SRC_OUTPUT_RESULT_VISUALIZATION_VISUALIZATION_H_

#include "output/result/resultstore.h"

namespace espreso {

class Visualization: public ResultStoreBase {

public:
	static bool storeStep();

	Visualization(const Mesh &mesh);
};

}



#endif /* SRC_OUTPUT_RESULT_VISUALIZATION_VISUALIZATION_H_ */
