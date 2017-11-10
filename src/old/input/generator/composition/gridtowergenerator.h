
#ifndef SRC_INPUT_GENERATOR_COMPOSITION_GRIDTOWERGENERATOR_H_
#define SRC_INPUT_GENERATOR_COMPOSITION_GRIDTOWERGENERATOR_H_

#include "gridgenerator.h"

namespace espreso {

struct GridTowerGeneratorConfiguration;

namespace input {

class GridTowerGenerator: public OldLoader {

public:
	GridTowerGenerator(const GridTowerGeneratorConfiguration &configuration, OldMesh &mesh, size_t index, size_t size);
	virtual ~GridTowerGenerator();

	static void load(const GridTowerGeneratorConfiguration &configuration, OldMesh &mesh, size_t index, size_t size);

	virtual void points(Coordinates &coordinates);
	virtual void elements(std::vector<size_t> &bodies, std::vector<OldElement*> &elements, std::vector<OldElement*> &faces, std::vector<OldElement*> &edges);
	virtual void neighbours(std::vector<OldElement*> &nodes, std::vector<int> &neighbours, const std::vector<OldElement*> &faces, const std::vector<OldElement*> &edges);
	virtual void regions(
			std::vector<Evaluator*> &evaluators,
			std::vector<Region*> &regions,
			std::vector<OldElement*> &elements,
			std::vector<OldElement*> &faces,
			std::vector<OldElement*> &edges,
			std::vector<OldElement*> &nodes);

	virtual bool partitiate(const std::vector<OldElement*> &nodes, std::vector<eslocal> &partsPtrs, std::vector<std::vector<OldElement*> > &fixPoints, std::vector<OldElement*> &corners);

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
