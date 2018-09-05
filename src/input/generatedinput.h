
#ifndef SRC_INPUT_GENERATEDINPUT_H_
#define SRC_INPUT_GENERATEDINPUT_H_

#include "input.h"

namespace espreso {

class GeneratedInput: public Input {

public:
	static void buildMesh(PlainMeshData &meshData, Mesh &mesh, bool needSynchronization=false);

protected:
	GeneratedInput(PlainMeshData &dMesh, Mesh &mesh, bool needSynchronization);

	void removeDanglingNodes();
	void synchronizeGlobalIndices();
};

}


#endif /* SRC_INPUT_GENERATEDINPUT_H_ */
