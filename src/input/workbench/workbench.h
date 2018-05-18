
#ifndef SRC_INPUT_WORKBENCH_WORKBENCH_H_
#define SRC_INPUT_WORKBENCH_WORKBENCH_H_

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

struct Point;
struct PlainMeshData;
class ECFRoot;
class Mesh;

class WorkbenchLoader {

public:
	static void load(const ECFRoot &configuration, Mesh &mesh);

protected:
	WorkbenchLoader(const ECFRoot &configuration, Mesh &mesh);

	void readData();
	void prepareData();
	void parseData(PlainMeshData &dMesh);

	const ECFRoot &_configuration;
	Mesh &_mesh;

	std::vector<NBlock> _NBlocks;
	std::vector<EBlock> _EBlocks;
//	std::vector<EBlock> _BBlocks;
	std::vector<CMBlock> _CMBlocks;
	std::vector<ET> _ET;
	std::vector<ESel> _ESel;
	std::vector<NSel> _NSel;
	std::vector<CM> _CM;
	std::vector<BlockEnd> _blockEnds;

	const char *_begin, *_current, *_end;
	std::vector<char> _data;
	std::vector<size_t> _dataOffset;
};

}



#endif /* SRC_INPUT_WORKBENCH_WORKBENCH_H_ */
