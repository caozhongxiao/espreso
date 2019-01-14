
#ifndef SRC_INPUT_OPENFOAM_OPENFOAM_H_
#define SRC_INPUT_OPENFOAM_OPENFOAM_H_

#include "input/plaindata.h"
#include "input/mpiloader/mpiloader.h"

#include <cstddef>
#include <string>
#include <vector>

namespace espreso {

class InputConfiguration;
class Mesh;
struct OpenFOAMSet;

struct PlainOpenFOAMData: public PlainMeshData {
	esint nelements;
	std::vector<esint> fIDs, fsize, fnodes, owner, neighbour;
};

class OpenFOAMLoader {

public:
	static void load(const InputConfiguration &configuration, Mesh &mesh);

protected:
	OpenFOAMLoader(const InputConfiguration &configuration, Mesh &mesh);

	void distributedReader(const std::string &file, ParallelFile &pfile, bool isMandatory);

	void readData();
	void parseData(PlainOpenFOAMData &mesh);

	void buildFaces(PlainOpenFOAMData &mesh);

	void collectFaces(PlainOpenFOAMData &mesh);
	void buildElements(PlainOpenFOAMData &mesh);

	const InputConfiguration &_configuration;
	MPISubset _loaders;

	ParallelFile _points, _faces, _owner, _neighbour, _boundary;
	ParallelFile _pointZones, _faceZones, _cellZones;

	std::vector<OpenFOAMSet> _sets;
	std::vector<esint> _fdist, _edist;
};

}



#endif /* SRC_INPUT_OPENFOAM_OPENFOAM_H_ */
