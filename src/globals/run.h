
#ifndef SRC_GLOBALS_RUN_H_
#define SRC_GLOBALS_RUN_H_

namespace espreso {

class ECFRoot;
class Mesh;
class DataHolder;

struct run {

	static ECFRoot *ecf;
	static Mesh *mesh;
	static DataHolder *data;

	static void storeMesh();
	static void storeSolution();
};

}

#endif /* SRC_GLOBALS_RUN_H_ */
