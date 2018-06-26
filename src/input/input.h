
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

struct PlainElement {
	eslocal id;
	int etype;
	int body;
	int material;
};

struct PlainMeshData {
	std::vector<eslocal> nIDs;
	std::vector<Point> coordinates;

	std::vector<eslocal> esize, enodes;
	std::vector<PlainElement> edata;

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

	void fillElements();

	const ECFRoot &_configuration;
	PlainMeshData &_meshData;
	Mesh &_mesh;

	std::vector<size_t> _nDistribution, _eDistribution;
};

}



#endif /* SRC_INPUT_INPUT_H_ */
