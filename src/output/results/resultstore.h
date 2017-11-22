
#ifndef SRC_OUTPUT_RESULTS_RESULTSTORE_H_
#define SRC_OUTPUT_RESULTS_RESULTSTORE_H_

namespace espreso {

struct OutputConfiguration;

class Store {

public:
	virtual void storePreprocessedData() =0;

	virtual void updateMesh() =0;
	virtual void updateSolution() =0;

	virtual ~Store() {};
};

class ResultStore: public Store {

public:
	ResultStore(const OutputConfiguration &configuration): _configuration(configuration) {};

	virtual void storePreprocessedData() {}

	virtual void updateMesh() {}
	virtual void updateSolution() {}

protected:
	const OutputConfiguration &_configuration;
};

}



#endif /* SRC_OUTPUT_RESULTS_RESULTSTORE_H_ */
