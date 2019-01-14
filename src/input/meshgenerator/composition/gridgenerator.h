
#ifndef SRC_INPUT_GENERATOR_COMPOSITION_GRIDGENERATOR_H_
#define SRC_INPUT_GENERATOR_COMPOSITION_GRIDGENERATOR_H_

#include "input/meshgenerator/primitives/gridsettings.h"

namespace espreso {

class Mesh;
struct PlainMeshData;
struct GridGeneratorConfiguration;
struct SphereGeneratorConfiguration;

class GridGenerator {

	friend class GridTowerGenerator;
	friend class SphereGenerator;

public:
	static void generate(const GridGeneratorConfiguration &configuration, Mesh &mesh);

protected:
	GridGenerator(const GridGeneratorConfiguration &configuration);
	GridGenerator(const SphereGeneratorConfiguration &configuration);
	virtual ~GridGenerator() {}

	virtual void init();
	virtual void nodes(PlainMeshData &mesh);
	virtual void elements(PlainMeshData &mesh);
	virtual void neighbors(PlainMeshData &mesh);
	virtual void regions(const GridGeneratorConfiguration &configuration, PlainMeshData &mesh);

	GridSettings _settings;
	BlockGenerator _block;
	Triple<size_t> _clusterOffset;
	std::vector<int> _clusterIndices;
};

}


#endif /* SRC_INPUT_GENERATOR_COMPOSITION_GRIDGENERATOR_H_ */
