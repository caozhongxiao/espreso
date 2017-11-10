
#ifndef SRC_INPUT_LOADER_H_
#define SRC_INPUT_LOADER_H_

namespace espreso {

class ECFConfiguration;
class Mesh;

class Loader {

public:
	static void load(const ECFConfiguration &configuration, Mesh &mesh, int MPIrank, int MPIsize);
};

}

#endif /* SRC_INPUT_LOADER_H_ */
