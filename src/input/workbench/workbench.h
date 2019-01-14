
#ifndef SRC_INPUT_WORKBENCH_WORKBENCH_H_
#define SRC_INPUT_WORKBENCH_WORKBENCH_H_

#include "input/plaindata.h"
#include "input/mpiloader/mpiloader.h"

#include <cstddef>
#include <vector>

namespace espreso {

class InputConfiguration;
class Mesh;

struct NBlock;
struct EBlock;
struct CMBlock;
struct ET;
struct ESel;
struct NSel;
struct CM;
struct BlockEnd;

struct PlainWorkbenchData: public PlainMeshData {
	std::vector<int> et;
};

class WorkbenchLoader {

public:
	static void load(const InputConfiguration &configuration, Mesh &mesh);

protected:
	WorkbenchLoader(const InputConfiguration &configuration, Mesh &mesh);

	void readData();
	void prepareData();
	void parseData(PlainWorkbenchData &dMesh);

	const InputConfiguration &_configuration;

	std::vector<NBlock> _NBlocks;
	std::vector<EBlock> _EBlocks;
	std::vector<CMBlock> _CMBlocks;
	std::vector<ET> _ET;
	std::vector<ESel> _ESel;
	std::vector<NSel> _NSel;
	std::vector<CM> _CM;
	std::vector<BlockEnd> _blockEnds;

	ParallelFile _pfile;
};

}



#endif /* SRC_INPUT_WORKBENCH_WORKBENCH_H_ */
