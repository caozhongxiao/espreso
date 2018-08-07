
#ifndef SRC_OUTPUT_RESULT_VISUALIZATION_INSITUVISUALIZATION_H_
#define SRC_OUTPUT_RESULT_VISUALIZATION_INSITUVISUALIZATION_H_

#include "visualization.h"

namespace espreso {

class InSituWrapper;

struct InSituVisualization: public Visualization {

	InSituVisualization(const Mesh &mesh, const OutputConfiguration &configuration);
	~InSituVisualization();

	virtual bool isCollected() { return false; }
	virtual bool isDistributed() { return true; }

	void updateMesh();
	void updateSolution(const Step &step);

protected:
	InSituWrapper *_inSitu;
};

}



#endif /* SRC_OUTPUT_RESULT_VISUALIZATION_INSITUVISUALIZATION_H_ */
