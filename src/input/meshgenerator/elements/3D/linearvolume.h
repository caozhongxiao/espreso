
#ifndef SRC_INPUT_MESHGENERATOR_ELEMENTS_3D_LINEARVOLUME_H_
#define SRC_INPUT_MESHGENERATOR_ELEMENTS_3D_LINEARVOLUME_H_

#include "../element.h"

namespace espreso {

struct LinearVolumeGenerator: public ElementGenerator {

	LinearVolumeGenerator();

	void pushNodes(std::vector<esint> &nodes, const std::vector<esint> &indices, CubeEdge edge) const;
	void pushTriangleNodes(std::vector<esint> &nodes, const std::vector<esint> &indices, CubeFace face) const;
	void pushSquareNodes(std::vector<esint> &nodes, const std::vector<esint> &indices, CubeFace face) const;

	void pushEdge(std::vector<esint> &elements, std::vector<esint> &esize, std::vector<int> &etype, const std::vector<esint> &indices, CubeEdge edge) const;
};

}


#endif /* SRC_INPUT_MESHGENERATOR_ELEMENTS_3D_LINEARVOLUME_H_ */
