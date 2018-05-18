
#ifndef SRC_INPUT_SORTEDINPUT_H_
#define SRC_INPUT_SORTEDINPUT_H_

#include "input.h"

namespace espreso {

class SortedInput: public Input {

public:
	static void load(const ECFRoot &configuration, Mesh &mesh, int MPIrank, int MPIsize);
	static void buildMesh(const ECFRoot &configuration, PlainMeshData &meshData, Mesh &mesh);

protected:
	SortedInput(const ECFRoot &configuration, PlainMeshData &dMesh, Mesh &mesh);

	void checkERegions();

	void fillElements();
	void fillCoordinates();
	void addNodeRegions();
	void addBoundaryRegions();
	void addElementRegions();

	std::vector<int> _targetRanks;
	std::vector<std::vector<eslocal> > _rankNodeMap;
};

}

#endif /* SRC_INPUT_SORTEDINPUT_H_ */
