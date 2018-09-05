
#ifndef SRC_INPUT_MESHGENERATOR_COMPOSITION_SPHEREGENERATOR_H_
#define SRC_INPUT_MESHGENERATOR_COMPOSITION_SPHEREGENERATOR_H_

#include "gridgenerator.h"

namespace espreso {

class Mesh;
struct PlainMeshData;
struct SphereGeneratorConfiguration;

class SphereGenerator: public GridGenerator {

public:
	static void generate(const SphereGeneratorConfiguration &configuration, Mesh &mesh);

protected:
	SphereGenerator(const SphereGeneratorConfiguration &configuration);
	virtual ~SphereGenerator() {}

	virtual void nodes(PlainMeshData &mesh);
	virtual void neighbors(PlainMeshData &mesh);
	virtual void regions(const SphereGeneratorConfiguration &configuration, PlainMeshData &mesh);

	enum class SIDE : int {
		UP = 0,
		FRONT = 1,
		DOWN = 2,
		BACK = 3,
		LEFT = 4,
		RIGHT = 5
	};

//	SphereSettings _settings;

	Triple<size_t> _clusterOffset;
	Triple<size_t> _subnodes;

	SIDE _side;
	size_t _row, _col, _layer;
};

}



#endif /* SRC_INPUT_MESHGENERATOR_COMPOSITION_SPHEREGENERATOR_H_ */
