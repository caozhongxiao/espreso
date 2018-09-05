
#ifndef SRC_INPUT_MESHGENERATOR_ELEMENTS_3D_PYRAMID5_H_
#define SRC_INPUT_MESHGENERATOR_ELEMENTS_3D_PYRAMID5_H_

#include "../element.h"

namespace espreso {

struct Pyramid5Generator: public ElementGenerator {

	Pyramid5Generator();

	void pushElements(std::vector<eslocal> &elements, const std::vector<eslocal> &indices) const;
	void pushNodes(std::vector<eslocal> &nodes, const std::vector<eslocal> &indices, CubeEdge edge) const;
	void pushNodes(std::vector<eslocal> &nodes, const std::vector<eslocal> &indices, CubeFace face) const;
	void pushEdge(std::vector<eslocal> &elements, std::vector<eslocal> &esize, std::vector<eslocal> &etype, const std::vector<eslocal> &indices, CubeEdge edge) const;
	void pushFace(std::vector<eslocal> &elements, std::vector<eslocal> &esize, std::vector<eslocal> &etype, const std::vector<eslocal> &indices, CubeFace face) const;
};

}

#endif /* SRC_INPUT_MESHGENERATOR_ELEMENTS_3D_PYRAMID5_H_ */
