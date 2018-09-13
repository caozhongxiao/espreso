
#ifndef SRC_INPUT_MESHGENERATOR_ELEMENTS_3D_HEXAHEDRON20_H_
#define SRC_INPUT_MESHGENERATOR_ELEMENTS_3D_HEXAHEDRON20_H_

#include "quadraticvolume.h"

namespace espreso {

struct Hexahedron20Generator: public QuadraticVolumeGenerator {

	Hexahedron20Generator();

	void pushElements(std::vector<eslocal> &elements, const std::vector<eslocal> &indices) const;
	void pushFace(std::vector<eslocal> &elements, std::vector<eslocal> &esize, std::vector<int> &etype, const std::vector<eslocal> &indices, CubeFace face) const;

	void pushNodes(std::vector<eslocal> &nodes, const std::vector<eslocal> &indices, CubeFace face) const
	{
		pushSquareNodes(nodes, indices, face);
	}
};

}


#endif /* SRC_INPUT_MESHGENERATOR_ELEMENTS_3D_HEXAHEDRON20_H_ */
