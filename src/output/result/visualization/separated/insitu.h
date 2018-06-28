
#ifndef SRC_OUTPUT_RESULT_VISUALIZATION_SEPARATED_INSITU_H_
#define SRC_OUTPUT_RESULT_VISUALIZATION_SEPARATED_INSITU_H_

#include "../visualization.h"

namespace espreso {

class InSituWrapper;

struct InSitu: public Visualization {

	InSitu(const Mesh &mesh, const OutputConfiguration &configuration);
	~InSitu();

	virtual bool isCollected() { return false; }
	virtual bool isSeparated() { return true; }

	void updateMesh();
	void updateSolution(const Step &step);

protected:
	InSituWrapper *_inSitu;
};

}



#endif /* SRC_OUTPUT_RESULT_VISUALIZATION_SEPARATED_INSITU_H_ */
