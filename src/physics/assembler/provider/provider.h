
#ifndef SRC_PHYSICS_ASSEMBLER_PROVIDER_PROVIDER_H_
#define SRC_PHYSICS_ASSEMBLER_PROVIDER_PROVIDER_H_

namespace espreso {

enum class MatrixType;
struct DataHolder;
struct LoadStepConfiguration;

class Provider {

public:
	Provider(DataHolder *data, LoadStepConfiguration &configuration);
	virtual ~Provider() {}

	virtual MatrixType getMatrixType() const =0;

	virtual bool needMatrixVectorProduct();
	virtual bool needOriginalStiffnessMatrices();
	virtual bool needOriginalRHS();
	virtual bool needReactionForces();
	virtual double& solutionPrecision() =0;

protected:
	DataHolder *_data;
	LoadStepConfiguration &_configuration;
};


}



#endif /* SRC_PHYSICS_ASSEMBLER_PROVIDER_PROVIDER_H_ */
