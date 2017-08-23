
#ifndef SRC_INPUT_GENERATOR_COMPOSITION_GRIDTOWERGENERATOR_H_
#define SRC_INPUT_GENERATOR_COMPOSITION_GRIDTOWERGENERATOR_H_

#include "gridgenerator.h"

namespace espreso {

struct GridTowerGeneratorConfiguration;

namespace input {

class GridTowerGenerator: public Loader {

public:
	GridTowerGenerator(const GridTowerGeneratorConfiguration &configuration, Mesh &mesh, size_t index, size_t size);
	virtual ~GridTowerGenerator();

	static void load(const GridTowerGeneratorConfiguration &configuration, Mesh &mesh, size_t index, size_t size);

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
	const GridTowerGeneratorConfiguration &_configuration;
	GridGenerator* _gridGenerator;

	size_t _clusterIndexBegin;
	size_t _clusterIndexEnd;
	size_t _gridIndex;
	size_t _gridPointsIDOffset;
	Triple<esglobal> _gridPointsOffset;

	size_t _index;
	size_t _size;
};


}
}



#endif /* SRC_INPUT_GENERATOR_COMPOSITION_GRIDTOWERGENERATOR_H_ */
