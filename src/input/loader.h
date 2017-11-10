
#ifndef SRC_INPUT_LOADER_H_
#define SRC_INPUT_LOADER_H_

namespace espreso {

class ECFConfiguration;
class NewMesh;

class Loader {

public:
	static void load(const ECFConfiguration &configuration, NewMesh &mesh, int MPIrank, int MPIsize);
};

}

#endif /* SRC_INPUT_LOADER_H_ */
