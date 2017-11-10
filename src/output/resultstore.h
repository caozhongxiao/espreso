
#ifndef SRC_OUTPUT_RESULTSTORE_H_
#define SRC_OUTPUT_RESULTSTORE_H_

namespace espreso {

class ResultStore {
public:
	virtual void updateMesh() {};
	virtual void updateSolution() {};

	ResultStore() {};
	virtual ~ResultStore() {};
};

}



#endif /* SRC_OUTPUT_RESULTSTORE_H_ */
