
#ifndef SRC_OUTPUT_RESULT_MONITORS_MONITORING_H_
#define SRC_OUTPUT_RESULT_MONITORS_MONITORING_H_

#include "mpi.h"

#include "../resultstore.h"
#include "../../../mesh/store/statisticsstore.h"

#include <utility>
#include <vector>
#include <fstream>

namespace espreso {

struct ElementsRegionStore;
struct BoundaryRegionStore;
struct NodeData;
struct ElementData;

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
	static bool storeStep(const OutputConfiguration &configuration, const Step &step);

	virtual bool isCollected() { return true; }
	virtual bool isSeparated() { return false; }

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

	std::vector<std::pair<NodeData*, const ElementsRegionStore*> > _nedata;
	std::vector<std::pair<NodeData*, const BoundaryRegionStore*> > _nbdata;
	std::vector<std::pair<ElementData*, const ElementsRegionStore*> > _edata;

	std::vector<Statistics> _data;
};

}



#endif /* SRC_OUTPUT_RESULT_MONITORS_MONITORING_H_ */
