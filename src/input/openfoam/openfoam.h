
#ifndef SRC_INPUT_OPENFOAM_OPENFOAM_H_
#define SRC_INPUT_OPENFOAM_OPENFOAM_H_

#include "../plaindata.h"
#include "../mpiloader/mpiloader.h"

#include <cstddef>
#include <vector>

namespace espreso {

class InputConfiguration;
class Mesh;
struct OpenFOAMBoundaryData;

struct PlainOpenFOAMData: public PlainMeshData {
	std::vector<eslocal> fID, fsize, fnodes, owner, neighbour;
};

class OpenFOAMLoader {

public:
	static void load(const InputConfiguration &configuration, Mesh &mesh);

protected:
	OpenFOAMLoader(const InputConfiguration &configuration, Mesh &mesh);

	void readData();
	void parseData(PlainOpenFOAMData &mesh);
	void collectFaces(PlainOpenFOAMData &mesh);
	void buildElements(PlainOpenFOAMData &mesh);
	void buildFaces(PlainOpenFOAMData &mesh);

	const InputConfiguration &_configuration;

	ParallelFile _points, _faces, _owner, _neighbour, _boundary;

	std::vector<OpenFOAMBoundaryData> _boundaryData;
	std::vector<eslocal> _fdist, _edist;
};

}



#endif /* SRC_INPUT_OPENFOAM_OPENFOAM_H_ */
