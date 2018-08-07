
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
	static void load(const ECFRoot &configuration, Mesh &mesh);

protected:
	Input(const ECFRoot &configuration, PlainMeshData &meshData, Mesh &mesh)
	: _configuration(configuration), _meshData(meshData), _mesh(mesh), _eregsize(1), _nregsize(1) {}

	void balance();

	std::vector<size_t> getDistribution(const std::vector<eslocal> &IDs, const std::vector<eslocal> &permutation);

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
	void sortElements(const std::vector<eslocal> &permutation);

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
