
#ifndef SRC_PHYSICS_ASSEMBLER_PROVIDER_PROVIDER_H_
#define SRC_PHYSICS_ASSEMBLER_PROVIDER_PROVIDER_H_

namespace espreso {

enum class MatrixType;
struct LoadStepConfiguration;

class Provider {

public:
	Provider(LoadStepConfiguration &configuration);
	virtual ~Provider() {}

	virtual MatrixType getMatrixType() const =0;

	virtual bool needOriginalStiffnessMatrices();
	virtual double& solutionPrecision() =0;

protected:
	LoadStepConfiguration &_configuration;
};


}



#endif /* SRC_PHYSICS_ASSEMBLER_PROVIDER_PROVIDER_H_ */
