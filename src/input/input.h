
#ifndef SRC_INPUT_INPUT_H_
#define SRC_INPUT_INPUT_H_

#include "plaindata.h"

#include <cstddef>
#include <string>
#include <vector>

namespace espreso {

class ECFRoot;
class Mesh;

class Input {

public:
	static bool load(const ECFRoot &configuration, Mesh &mesh);

protected:
	Input(PlainMeshData &meshData, Mesh &mesh)
	: _meshData(meshData), _mesh(mesh), _eregsize(1), _nregsize(1) {}

	void balance();

	std::vector<esint> getDistribution(const std::vector<esint> &IDs, const std::vector<esint> &permutation);

	void balanceNodes();
	void balancePermutedNodes();
	void balanceElements();
	void balancePermutedElements();

	void assignRegions(
			std::map<std::string, std::vector<esint> > &regions, std::vector<esint> &IDs,
			std::vector<esint> &distribution,
			size_t &rsize, std::vector<esint> &rbits);
	void fillRegions(std::map<std::string, std::vector<esint> > &regions, size_t &rsize, std::vector<esint> &rbits);

	void sortNodes(bool withElementNodes = false);
	void sortElements();
	void sortElements(const std::vector<esint> &permutation);

	void fillNodes();
	void fillElements();
	void fillNeighbors();

	void fillNodeRegions();
	void fillBoundaryRegions();
	void fillElementRegions();

	void reindexRegions(std::map<std::string, std::vector<esint> > &regions, std::vector<esint> &IDs);
	void reindexElementNodes();
	void reindexBoundaryNodes();

	PlainMeshData &_meshData;
	Mesh &_mesh;

	std::vector<esint> _nDistribution, _eDistribution, _etypeDistribution;

	size_t _eregsize, _nregsize;
	std::vector<esint> _eregions, _nregions;
};

}



#endif /* SRC_INPUT_INPUT_H_ */
