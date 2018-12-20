
#ifndef SRC_PHYSICS_ASSEMBLER_PROVIDER_PROVIDER_H_
#define SRC_PHYSICS_ASSEMBLER_PROVIDER_PROVIDER_H_

namespace espreso {

struct LoadStepConfiguration;

class Provider {

public:
	Provider(LoadStepConfiguration &configuration);
	virtual ~Provider() {}

	virtual bool needOriginalStiffnessMatrices();

protected:
	LoadStepConfiguration &_configuration;
};


}



#endif /* SRC_PHYSICS_ASSEMBLER_PROVIDER_PROVIDER_H_ */
