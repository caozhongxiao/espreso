
#ifndef SRC_INPUT_MESHGENERATOR_ELEMENTS_3D_TETRAHEDRON10_H_
#define SRC_INPUT_MESHGENERATOR_ELEMENTS_3D_TETRAHEDRON10_H_

#include "quadraticvolume.h"

namespace espreso {

struct Tetrahedron10Generator: public QuadraticVolumeGenerator {

	Tetrahedron10Generator();

	void pushElements(std::vector<eslocal> &elements, const std::vector<eslocal> &indices) const;
	void pushFace(std::vector<eslocal> &elements, std::vector<eslocal> &esize, std::vector<int> &etype, const std::vector<eslocal> &indices, CubeFace face) const;

	void pushNodes(std::vector<eslocal> &nodes, const std::vector<eslocal> &indices, CubeFace face) const
	{
		pushTriangleNodes(nodes, indices, face);
	}
};

}


#endif /* SRC_INPUT_MESHGENERATOR_ELEMENTS_3D_TETRAHEDRON10_H_ */
