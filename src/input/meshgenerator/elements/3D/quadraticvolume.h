
#ifndef SRC_INPUT_MESHGENERATOR_ELEMENTS_3D_QUADRATICVOLUME_H_
#define SRC_INPUT_MESHGENERATOR_ELEMENTS_3D_QUADRATICVOLUME_H_

#include "../element.h"

namespace espreso {

struct QuadraticVolumeGenerator: public ElementGenerator {

	QuadraticVolumeGenerator();

	void pushNodes(std::vector<eslocal> &nodes, const std::vector<eslocal> &indices, CubeEdge edge) const;
	void pushTriangleNodes(std::vector<eslocal> &nodes, const std::vector<eslocal> &indices, CubeFace face) const;
	void pushSquareNodes(std::vector<eslocal> &nodes, const std::vector<eslocal> &indices, CubeFace face) const;

	void pushEdge(std::vector<eslocal> &elements, std::vector<eslocal> &esize, std::vector<eslocal> &etype, const std::vector<eslocal> &indices, CubeEdge edge) const;
};

}

#endif /* SRC_INPUT_MESHGENERATOR_ELEMENTS_3D_QUADRATICVOLUME_H_ */
