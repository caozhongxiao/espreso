
#ifndef SRC_INPUT_PLAINDATA_H_
#define SRC_INPUT_PLAINDATA_H_

#include "../basis/containers/point.h"

#include <vector>
#include <map>

namespace espreso {

// Plain mesh data for building the mesh in ESPRESO
// parameters with prefix '_' are optional -> loaders fill them if neccessary

struct PlainMeshData {
	std::vector<eslocal> nIDs;       // nodes IDs [arbitrary numbers]
	std::vector<Point> coordinates;  // nodes coordinates

	std::vector<eslocal> _nrankdist; // nodes ranks distribution [0, 2, 5, ...] n0 is on 2 processes
	std::vector<int> _nranks;        // nodes ranks              [0, 1, ...]    n0 is on processes 0 and 1


	std::vector<eslocal> eIDs;       // elements IDs [arbitrary numbers]
	std::vector<eslocal> esize;      // the number of nodes for a given element [4, 4, 4]         e0 has 4 nodes
	std::vector<eslocal> enodes;     // elements nodes                          [0, 1, 5, 4, ...] e0 has nodes 0, 1, 5, 4
	std::vector<int> etype;          // elements types [from Element::CODE]
	std::vector<int> body;           // elements bodies
	std::vector<int> material;       // elements materials

	std::vector<eslocal> _edist;     // elements nodes distribution [0, 4, 8, ...]


	std::map<std::string, std::vector<eslocal> > eregions; // elements regions <name, list of IDs>
	std::map<std::string, std::vector<eslocal> > nregions; // nodes regions <name, list of IDs>
};

}



#endif /* SRC_INPUT_PLAINDATA_H_ */
