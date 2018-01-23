
#ifndef SRC_INPUT_WORKBENCH_WORKBENCH_H_
#define SRC_INPUT_WORKBENCH_WORKBENCH_H_

#include "parser/nblock.h"
#include "parser/eblock.h"
#include "parser/cmblock.h"
#include "parser/blockend.h"

#include <vector>

namespace espreso {

struct Point;
struct DistributedMesh;
class ECFConfiguration;
class Mesh;

class WorkbenchLoader {

public:
	static void load(const ECFConfiguration &configuration, Mesh &mesh);

protected:
	WorkbenchLoader(const ECFConfiguration &configuration, Mesh &mesh);

	void readData();
	void prepareData();
	void parseData(DistributedMesh &dMesh);

	const ECFConfiguration &_configuration;
	Mesh &_mesh;

	std::vector<NBlock> _NBlocks;
	std::vector<EBlock> _EBlocks;
	std::vector<EBlock> _BBlocks;
	std::vector<CMBlock> _CMBlocks;
	std::vector<BlockEnd> _blockEnds;

	const char *_begin, *_current, *_end;
	std::vector<char> _data;
	std::vector<eslocal> _dataOffset;
};

}



#endif /* SRC_INPUT_WORKBENCH_WORKBENCH_H_ */
