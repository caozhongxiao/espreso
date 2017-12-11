
#ifndef SRC_OUTPUT_RESULT_MONITORS_MONITORING_H_
#define SRC_OUTPUT_RESULT_MONITORS_MONITORING_H_

#include "../resultstore.h"

#include <vector>
#include <fstream>

namespace espreso {

struct Monitor {
	eslocal printSize;
//	Region* region;
//	std::vector<Property> properties;

	Monitor();
};

class Monitoring: public ResultStoreBase {

public:
	void updateMesh() {}
	void updateSolution(const Step &step);

	Monitoring(const Mesh &mesh, const OutputConfiguration &configuration);
	static char delimiter;

protected:
	const Mesh &_mesh;
	std::ofstream _os;

	std::vector<Monitor> _monitors;
};

}



#endif /* SRC_OUTPUT_RESULT_MONITORS_MONITORING_H_ */
