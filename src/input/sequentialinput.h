
#ifndef SRC_INPUT_SEQUENTIALINPUT_H_
#define SRC_INPUT_SEQUENTIALINPUT_H_

#include "input.h"

namespace espreso {

class SequentialInput: public Input {

public:
	static void buildMesh(const ECFRoot &configuration, PlainMeshData &meshData, Mesh &mesh);

protected:
	SequentialInput(const ECFRoot &configuration, PlainMeshData &dMesh, Mesh &mesh);
};

}


#endif /* SRC_INPUT_SEQUENTIALINPUT_H_ */
