
#ifndef SRC_OUTPUT_SOLUTION_SOLUTIONSTORE_H_
#define SRC_OUTPUT_SOLUTION_SOLUTIONSTORE_H_

namespace espreso {

struct OutputConfiguration;
class Mesh;
struct Step;

class Store {

public:
	virtual void storePreprocessedData() =0;

	virtual bool storeSolution(const Step &step) =0;
	virtual bool storeStatistics(const Step &step) =0;

	virtual void updateMesh() =0;
	virtual void updateSolution(const Step &step) =0;

	virtual ~Store() {};
};

class SolutionStore: public Store {

public:
	SolutionStore(const Mesh &mesh, const OutputConfiguration &configuration)
	: _mesh(mesh), _configuration(configuration) {};

	virtual void storePreprocessedData() { }

	virtual bool storeSolution(const Step &step) { return false; }
	virtual bool storeStatistics(const Step &step) { return false; }

	virtual void updateMesh() {}
	virtual void updateSolution(const Step &step) {}

protected:
	const Mesh &_mesh;
	const OutputConfiguration &_configuration;
};

}



#endif /* SRC_OUTPUT_SOLUTION_SOLUTIONSTORE_H_ */
