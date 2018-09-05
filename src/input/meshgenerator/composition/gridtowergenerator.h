
#ifndef SRC_INPUT_MESHGENERATOR_COMPOSITION_GRIDTOWERGENERATOR_H_
#define SRC_INPUT_MESHGENERATOR_COMPOSITION_GRIDTOWERGENERATOR_H_

#include "gridgenerator.h"

namespace espreso {

class Mesh;
struct PlainMeshData;
struct GridTowerGeneratorConfiguration;

class GridTowerGenerator: public GridGenerator {

public:
	static void generate(const GridTowerGeneratorConfiguration &configuration, Mesh &mesh);

	static size_t gridIndex(const GridTowerGeneratorConfiguration &configuration);

protected:
	GridTowerGenerator(const GridTowerGeneratorConfiguration &configuration);
	virtual ~GridTowerGenerator() {}

	virtual void init(const GridTowerGeneratorConfiguration &configuration);
	virtual void nodes(PlainMeshData &mesh);
	virtual void neighbors(const GridTowerGeneratorConfiguration &configuration, PlainMeshData &mesh);
	virtual void regions(const GridTowerGeneratorConfiguration &configuration, PlainMeshData &mesh);

	size_t _gridIndex;
	size_t _gridNodeOffset;
};

}


#endif /* SRC_INPUT_MESHGENERATOR_COMPOSITION_GRIDTOWERGENERATOR_H_ */
