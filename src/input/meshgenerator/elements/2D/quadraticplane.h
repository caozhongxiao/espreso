
#ifndef SRC_INPUT_MESHGENERATOR_ELEMENTS_2D_QUADRATICPLANE_H_
#define SRC_INPUT_MESHGENERATOR_ELEMENTS_2D_QUADRATICPLANE_H_

#include "../element.h"

namespace espreso {

struct QuadraticPlaneGenerator: public ElementGenerator {

	QuadraticPlaneGenerator();

	void pushNodes(std::vector<eslocal> &nodes, const std::vector<eslocal> &indices, CubeEdge edge) const;
	void pushNodes(std::vector<eslocal> &nodes, const std::vector<eslocal> &indices, CubeFace face) const;

	void pushEdge(std::vector<eslocal> &elements, std::vector<eslocal> &esize, std::vector<int> &etype, const std::vector<eslocal> &indices, CubeEdge edge) const;
	void pushFace(std::vector<eslocal> &elements, std::vector<eslocal> &esize, std::vector<int> &etype, const std::vector<eslocal> &indices, CubeFace face) const;
};

}



#endif /* SRC_INPUT_MESHGENERATOR_ELEMENTS_2D_QUADRATICPLANE_H_ */
