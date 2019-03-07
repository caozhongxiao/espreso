
#ifndef SRC_PHYSICS_ASSEMBLER_PROVIDER_PROVIDER_H_
#define SRC_PHYSICS_ASSEMBLER_PROVIDER_PROVIDER_H_

namespace espreso {

enum class MatrixType;
struct DataHolder;
struct LoadStepSolverConfiguration;

class Provider {

public:
	Provider(DataHolder *data, LoadStepSolverConfiguration &configuration);
	virtual ~Provider() {}

	virtual MatrixType getMatrixType() const =0;

	virtual bool needMatrixVectorProduct();
	virtual bool needOriginalStiffnessMatrices();
	virtual bool needSolverStiffnessMatrices();
	virtual bool needSolverRHS();
	virtual bool needReactionForces();
	virtual double& solutionPrecision() =0;

protected:
	DataHolder *_data;
	LoadStepSolverConfiguration &_configuration;
};


}



#endif /* SRC_PHYSICS_ASSEMBLER_PROVIDER_PROVIDER_H_ */
