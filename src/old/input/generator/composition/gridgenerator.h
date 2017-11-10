
#ifndef SRC_INPUT_GENERATOR_COMPOSITION_GRIDGENERATOR_H_
#define SRC_INPUT_GENERATOR_COMPOSITION_GRIDGENERATOR_H_

#include "../../loader.h"
#include "../primitives/triple.h"

#include "../../../basis/expression/expression.h"

namespace espreso {

enum class GENERATOR_ELEMENT_TYPE;
struct GridGeneratorConfiguration;

namespace input {

struct BlockGenerator;

struct GridSettings {

	GridSettings();
	GridSettings(const GridGeneratorConfiguration &configuration);

	GENERATOR_ELEMENT_TYPE etype;

	Triple<esglobal> start, end;
	Triple<size_t> blocks, clusters, domains, elements;
	Triple<Expression> projection, rotation;

	std::vector<bool> nonempty;
	bool uniformDecomposition;
};

class GridGenerator: public OldLoader {

	friend class GridTowerGenerator;

public:
	GridGenerator(const GridGeneratorConfiguration &configuration, OldMesh &mesh, size_t index, size_t size);
	virtual ~GridGenerator();

	static void load(const GridGeneratorConfiguration &configuration, OldMesh &mesh, size_t index, size_t size);

	virtual void points(Coordinates &coordinates) { points(coordinates, 0); }
	virtual void points(Coordinates &coordinates, size_t globalIdOffset);
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

	size_t pointCount() const;

	size_t bodyIndex() { return _body; }
	void bodyIndex(size_t index) { _body = index; }

protected:
	const GridGeneratorConfiguration &_configuration;
	GridSettings _settings;
	BlockGenerator* _block;
	Triple<size_t> _clusterOffset;
	Triple<size_t> _subnodes;

	size_t _index;
	size_t _size;
	size_t _body;

	std::vector<int> _cMap;
	int _clusterIndex;
};


}
}



#endif /* SRC_INPUT_GENERATOR_COMPOSITION_GRIDGENERATOR_H_ */
