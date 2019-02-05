
#ifndef SRC_INPUT_ABAQUS_ABAQUS_H_
#define SRC_INPUT_ABAQUS_ABAQUS_H_

#include "../plaindata.h"
#include "../mpiloader/mpiloader.h"

#include <cstddef>
#include <vector>

namespace espreso {

class InputConfiguration;
class Mesh;

struct NList;
struct EList;
struct CMBlock;
struct ET;
struct Eset;
struct Nset;
struct CM;
struct BlockFinish;
struct SSection;
struct Materials;
struct Elemat;

struct PlainAbaqusData: public PlainMeshData {
	std::vector<int> et;
};

class AbaqusLoader {

public:
	static void load(const InputConfiguration &configuration, Mesh &mesh);

protected:
	AbaqusLoader(const InputConfiguration &configuration, Mesh &mesh);

	void readData();
	void prepareData();
	void parseData(PlainAbaqusData &dMesh);
	const InputConfiguration &_configuration;

	std::vector<NList> _NLists;
	std::vector<EList> _ELists;
	std::vector<BlockFinish> _blockFinishs;
	std::vector<Eset> _Esets;
	std::vector<SSection> _SSections;
    std::vector<Materials> _Materials;
    std::vector<Elemat> _Elemats;
    std::vector<Nset> _Nsets;
	ParallelFile _pfile;
};

}



#endif /* SRC_INPUT_ABAQUS_ABAQUS_H_ */
