
#ifndef SRC_OUTPUT_RESULT_MONITORS_MONITORING_H_
#define SRC_OUTPUT_RESULT_MONITORS_MONITORING_H_

#include "mpi.h"

#include "../resultstore.h"
#include "../../../mesh/store/statisticsstore.h"

#include <vector>
#include <fstream>

namespace espreso {

struct ElementsRegionStore;
struct BoundaryRegionStore;

struct Monitor {
	std::string name;
	std::string property;
	std::string stats;
	eslocal printSize;
	double *data;

	Monitor(): name("---"), property("-"), printSize(5), data(NULL) {}
};

class Monitoring: public ResultStoreBase {

public:
	virtual bool isCollected() { return true; }

	void updateMesh();
	void updateSolution(const Step &step);

	Monitoring(const Mesh &mesh, const OutputConfiguration &configuration, bool async);
	~Monitoring();
	static char delimiter;

protected:
	const OutputConfiguration &_configuration;
	bool _async;
	MPI_Comm _communicator;
	std::ofstream _os;

	std::vector<Monitor> _monitors;

	std::vector<const ElementsRegionStore*> _eregions;
	std::vector<const BoundaryRegionStore*> _bregions;
	std::vector<Statistics> _data;
};

}



#endif /* SRC_OUTPUT_RESULT_MONITORS_MONITORING_H_ */
