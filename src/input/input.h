
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
	: _configuration(configuration), _meshData(meshData), _mesh(mesh), _eregsize(1), _nregsize(1) {}

	void balance();
	void balanceNodes();
	void balancePermutedNodes();
	void balanceElements();
	void balancePermutedElements();

	void assignRegions(
			std::map<std::string, std::vector<eslocal> > &regions, std::vector<eslocal> &IDs,
			std::vector<size_t> &distribution,
			size_t &rsize, std::vector<eslocal> &rbits);
	void fillRegions(std::map<std::string, std::vector<eslocal> > &regions, size_t &rsize, std::vector<eslocal> &rbits);

	void sortNodes();
	void sortElements();

	void fillNodes();
	void fillElements();

	void fillNodeRegions();
	void fillBoundaryRegions();
	void fillElementRegions();

	void reindexRegions(std::map<std::string, std::vector<eslocal> > &regions, std::vector<eslocal> &IDs);
	void reindexElementNodes();
	void reindexBoundaryNodes();

	const ECFRoot &_configuration;
	PlainMeshData &_meshData;
	Mesh &_mesh;

	std::vector<size_t> _nDistribution, _eDistribution, _etypeDistribution;

	size_t _eregsize, _nregsize;
	std::vector<eslocal> _eregions, _nregions;
};

}



#endif /* SRC_INPUT_INPUT_H_ */
