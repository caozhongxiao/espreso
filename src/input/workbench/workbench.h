
#ifndef SRC_INPUT_WORKBENCH_WORKBENCH_H_
#define SRC_INPUT_WORKBENCH_WORKBENCH_H_

#include "../plaindata.h"
#include "../mpiloader/mpiloader.h"

#include "parser/nblock.h"
#include "parser/eblock.h"
#include "parser/cmblock.h"
#include "parser/et.h"
#include "parser/esel.h"
#include "parser/nsel.h"
#include "parser/cm.h"
#include "parser/blockend.h"

#include <cstddef>
#include <vector>

namespace espreso {

class InputConfiguration;
class Mesh;

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

	ParallelFile pfile;
};

}



#endif /* SRC_INPUT_WORKBENCH_WORKBENCH_H_ */
