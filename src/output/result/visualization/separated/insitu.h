
#ifndef SRC_OUTPUT_RESULT_VISUALIZATION_SEPARATED_INSITU_H_
#define SRC_OUTPUT_RESULT_VISUALIZATION_SEPARATED_INSITU_H_

#include "output/result/visualization/visualization.h"

namespace espreso {

class Catalyst;

struct InSitu: public Visualization {

	InSitu(const Mesh &mesh, const OutputConfiguration &configuration);
	~InSitu();

	virtual bool isCollected() { return false; }
	virtual bool isSeparated() { return true; }

	void updateMesh();
	void updateSolution();

protected:
	Catalyst *_catalyst;
};

}



#endif /* SRC_OUTPUT_RESULT_VISUALIZATION_SEPARATED_INSITU_H_ */
