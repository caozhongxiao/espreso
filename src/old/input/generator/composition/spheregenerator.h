
#ifndef SRC_INPUT_GENERATOR_COMPOSITION_SPHEREGENERATOR_H_
#define SRC_INPUT_GENERATOR_COMPOSITION_SPHEREGENERATOR_H_

#include "../../loader.h"

#include "../primitives/triple.h"

namespace espreso {

enum class GENERATOR_ELEMENT_TYPE;
struct SphereGeneratorConfiguration;

namespace input {

struct BlockGenerator;

struct SphereSettings {

	SphereSettings();
	SphereSettings(const SphereGeneratorConfiguration &configuration);

	GENERATOR_ELEMENT_TYPE etype;

	double innerRadius, outerRadius;
	size_t clusters, layers;
	Triple<size_t> domains, elements;

	bool uniformDecomposition;
};

class SphereGenerator: public OldLoader {

public:
	SphereGenerator(const SphereGeneratorConfiguration &configuration, OldMesh &mesh, size_t index, size_t size);
	virtual ~SphereGenerator();

	static void load(const SphereGeneratorConfiguration &configuration, OldMesh &mesh, size_t index, size_t size);

	virtual void points(Coordinates &coordinates);
	virtual void elements(std::vector<size_t> &bodies, std::vector<Element*> &elements, std::vector<Element*> &faces, std::vector<Element*> &edges);
	virtual void neighbours(std::vector<Element*> &nodes, std::vector<int> &neighbours, const std::vector<Element*> &faces, const std::vector<Element*> &edges);
	virtual void regions(
			std::vector<Evaluator*> &evaluators,
			std::vector<Region*> &regions,
			std::vector<Element*> &elements,
			std::vector<Element*> &faces,
			std::vector<Element*> &edges,
			std::vector<Element*> &nodes);

	virtual bool partitiate(const std::vector<Element*> &nodes, std::vector<eslocal> &partsPtrs, std::vector<std::vector<Element*> > &fixPoints, std::vector<Element*> &corners);

protected:
	enum class SIDE : int {
		UP = 0,
		FRONT = 1,
		DOWN = 2,
		BACK = 3,
		LEFT = 4,
		RIGHT = 5
	};

	const SphereGeneratorConfiguration &_configuration;
	SphereSettings _settings;
	BlockGenerator* _block;
	Triple<size_t> _clusterOffset;
	Triple<size_t> _subnodes;

	size_t _index;
	size_t _size;

	SIDE _side;
	size_t _row, _col, _layer;
};


}
}


#endif /* SRC_INPUT_GENERATOR_COMPOSITION_SPHEREGENERATOR_H_ */
