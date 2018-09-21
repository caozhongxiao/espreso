
#ifndef SRC_INPUT_OPENFOAM_OPENFOAM_H_
#define SRC_INPUT_OPENFOAM_OPENFOAM_H_

#include "../plaindata.h"
#include "../mpiloader/mpiloader.h"

#include <cstddef>
#include <vector>

namespace espreso {

class InputConfiguration;
class Mesh;

struct PlainOpenFOAMData: public PlainMeshData {
	std::vector<eslocal> fsize, fnodes, owner, neighbour;
};

class OpenFOAMLoader {

public:
	static void load(const InputConfiguration &configuration, Mesh &mesh);

protected:
	OpenFOAMLoader(const InputConfiguration &configuration, Mesh &mesh);

	void readData();
	void parseData(PlainOpenFOAMData &mesh);
	void buildElements(PlainOpenFOAMData &mesh);

	const InputConfiguration &_configuration;

	ParallelFile points, faces, owner, neighbour, boundary;
};

}



#endif /* SRC_INPUT_OPENFOAM_OPENFOAM_H_ */
