
#ifndef SRC_OUTPUT_SOLUTION_VISUALIZATION_COLLECTEDVISUALIZATION_H_
#define SRC_OUTPUT_SOLUTION_VISUALIZATION_COLLECTEDVISUALIZATION_H_

#include "mpi.h"

namespace espreso {

struct CollectedVisualization {
	CollectedVisualization();
	~CollectedVisualization();

protected:
	MPI_Comm _storeCommunicator;
};

}



#endif /* SRC_OUTPUT_SOLUTION_VISUALIZATION_COLLECTEDVISUALIZATION_H_ */
