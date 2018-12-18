
#ifndef SRC_PHYSICS_ASSEMBLER_PROVIDER_PROVIDER_H_
#define SRC_PHYSICS_ASSEMBLER_PROVIDER_PROVIDER_H_

namespace espreso {

class Provider {

public:
	virtual ~Provider() {}

	virtual bool needOriginalStiffnessMatrices() =0;
};


}



#endif /* SRC_PHYSICS_ASSEMBLER_PROVIDER_PROVIDER_H_ */
