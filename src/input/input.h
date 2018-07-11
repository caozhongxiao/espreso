
#ifndef SRC_INPUT_INPUT_H_
#define SRC_INPUT_INPUT_H_

#include "../basis/containers/point.h"

#include <cstddef>
#include <string>
#include <vector>
#include <map>

namespace espreso {

class ECFRoot;
class Mesh;

struct PlainMeshData {
	std::vector<eslocal> nIDs, ndist, nranks;
	std::vector<Point> coordinates;

	std::vector<eslocal> eIDs, esize, enodes;
	std::vector<int> etype, body, material;

	std::map<std::string, std::vector<eslocal> > eregions;
	std::map<std::string, std::vector<eslocal> > nregions;
};

class Input {

public:
	static void load(const ECFRoot &configuration, Mesh &mesh);

protected:
	Input(const ECFRoot &configuration, PlainMeshData &meshData, Mesh &mesh)
	: _configuration(configuration), _meshData(meshData), _mesh(mesh) {}

	void balance();
	void balanceNodes();
	void balancePermutedNodes();
	void balanceElements();
	void balancePermutedElements();

	void assignRegions();

	void sortNodes();
	void sortNodesWithRegion();
	void sortElements();

	void fillNodes();
	void fillElements();

	void reindexRegions();
	void reindexNodeRegions();
	void reindexElementRegions();
	void reindexElementNodes();

	const ECFRoot &_configuration;
	PlainMeshData &_meshData;
	Mesh &_mesh;

	std::vector<size_t> _nDistribution, _eDistribution, _etypeDistribution;
};

}



#endif /* SRC_INPUT_INPUT_H_ */
