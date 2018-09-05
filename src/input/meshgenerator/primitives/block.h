
#ifndef SRC_INPUT_GENERATOR_PRIMITIVES_BLOCK_H_
#define SRC_INPUT_GENERATOR_PRIMITIVES_BLOCK_H_

#include "triple.h"

#include <vector>
#include <functional>

namespace espreso {

enum class Pattern {
	CHESSBOARD_WHITE,
	CHESSBOARD_BLACK
};

struct Point;
struct BlockGeneratorConfiguration;
struct BlockSettings;
struct ElementGenerator;
struct PlainMeshData;
struct BlockBorder;

class BlockGenerator {

public:
	BlockGenerator(const BlockGeneratorConfiguration &configuration, const BlockSettings &block);
	~BlockGenerator();

	const ElementGenerator* element() const { return _element; }

	void coordinates(PlainMeshData &mesh);
	void elements(PlainMeshData &mesh);
	void neighbors(const std::vector<int> &surroundings, PlainMeshData &mesh);

	void nodesRegion(const BlockBorder &border, std::vector<eslocal> &nodes);
	void edgesRegion(const BlockBorder &border, PlainMeshData &mesh, std::vector<eslocal> &elements);
	void facesRegion(const BlockBorder &border, PlainMeshData &mesh, std::vector<eslocal> &elements);
	void elementsRegion(const BlockBorder &border, std::vector<eslocal> &elements);
	void pattern(const Triple<size_t> &offset, const Triple<size_t> &size, std::vector<eslocal> &elements, Pattern pattern, size_t psize);

private:
	bool region(
			BlockBorder &intersection,
			Triple<size_t> &ebegin, Triple<size_t> &eend,
			Triple<size_t> &nbegin, Triple<size_t> &nend);

	void forEachElement(
			const Triple<size_t> &begin,
			const Triple<size_t> &end,
			std::function<void(std::vector<eslocal> &indices)> operation,
			std::function<void(Triple<size_t> &offset)> restriction);

	const BlockSettings &_block;
	const ElementGenerator *_element;
};

}

#endif /* SRC_INPUT_GENERATOR_PRIMITIVES_BLOCK_H_ */
