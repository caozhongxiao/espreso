
#ifndef SRC_INPUT_PLAINDATA_H_
#define SRC_INPUT_PLAINDATA_H_

#include "../basis/containers/point.h"

#include <vector>
#include <map>

namespace espreso {

struct PlainMeshData {
	std::vector<eslocal> nIDs, ndist;
	std::vector<int> nranks;
	std::vector<Point> coordinates;

	std::vector<eslocal> eIDs, esize, enodes;
	std::vector<int> etype, body, material;

	std::map<std::string, std::vector<eslocal> > eregions;
	std::map<std::string, std::vector<eslocal> > nregions;
};

}



#endif /* SRC_INPUT_PLAINDATA_H_ */
