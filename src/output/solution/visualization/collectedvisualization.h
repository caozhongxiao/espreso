
#ifndef SRC_OUTPUT_SOLUTION_VISUALIZATION_COLLECTEDVISUALIZATION_H_
#define SRC_OUTPUT_SOLUTION_VISUALIZATION_COLLECTEDVISUALIZATION_H_

#include "mpi.h"
#include "visualization.h"

namespace espreso {

struct CollectedVisualization: public Visualization {
	CollectedVisualization();
	~CollectedVisualization();

protected:
	MPI_Comm _storeCommunicator;
};

}



#endif /* SRC_OUTPUT_SOLUTION_VISUALIZATION_COLLECTEDVISUALIZATION_H_ */
