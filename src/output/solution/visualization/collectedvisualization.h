
#ifndef SRC_OUTPUT_SOLUTION_VISUALIZATION_COLLECTEDVISUALIZATION_H_
#define SRC_OUTPUT_SOLUTION_VISUALIZATION_COLLECTEDVISUALIZATION_H_

#include "mpi.h"
#include "visualization.h"

#include <vector>

namespace espreso {

struct CollectedVisualization: public Visualization {
	CollectedVisualization();
	~CollectedVisualization();

protected:
	void clearIntervals();
	void pushInterval(eslocal size);
	MPI_Datatype* commitIntervals();
	void storeIntervals(const std::string &name, const std::string &data, MPI_Datatype* datatype);

	MPI_Comm _storeCommunicator;
	std::vector<MPI_Datatype*> _datatypes;

	eslocal _loffset, _goffset, _lsize, _gsize;

	std::vector<MPI_Aint> _displacement;
	std::vector<int> _lenghts;
};

}



#endif /* SRC_OUTPUT_SOLUTION_VISUALIZATION_COLLECTEDVISUALIZATION_H_ */
