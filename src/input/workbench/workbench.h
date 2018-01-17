
#ifndef SRC_INPUT_WORKBENCH_WORKBENCH_H_
#define SRC_INPUT_WORKBENCH_WORKBENCH_H_

#include <string>
#include <vector>

namespace espreso {

class ECFConfiguration;
class Mesh;

struct DataInterval {
	eslocal header;
	eslocal indexsize;
	eslocal datacount, datasize;
	eslocal sIndex, eIndex;
	int sRank, eRank;

	DataInterval()
	: header(-1),
	  indexsize(-1),
	  datacount(-1), datasize(-1),
	  sIndex(-1), eIndex(-1), sRank(-1), eRank(-1) {}
};

class WorkbenchLoader {

public:
	static void load(const ECFConfiguration &configuration, Mesh &mesh);

protected:
	WorkbenchLoader(const ECFConfiguration &configuration, Mesh &mesh);

	void readData();
	void parseData();

	std::string getLine(eslocal index);

	int elementNodeCount(int etype);

	const ECFConfiguration &_configuration;
	Mesh &_mesh;

	char *_begin, *_current, *_end;
	std::vector<char> _data;

	std::vector<eslocal> _dataOffset;
	std::vector<DataInterval> _coordinates;
	std::vector<eslocal> _etypes;
	std::vector<DataInterval> _elements;
	std::vector<DataInterval> _cmblocks;
};

}



#endif /* SRC_INPUT_WORKBENCH_WORKBENCH_H_ */
