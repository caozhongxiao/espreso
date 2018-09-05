
#ifndef SRC_INPUT_MESHGENERATOR_ELEMENTS_3D_PRISMA15_H_
#define SRC_INPUT_MESHGENERATOR_ELEMENTS_3D_PRISMA15_H_

#include "quadraticvolume.h"

namespace espreso {

struct Prisma15Generator: public QuadraticVolumeGenerator {

	Prisma15Generator();

	void pushElements(std::vector<eslocal> &elements, const std::vector<eslocal> &indices) const;
	void pushFace(std::vector<eslocal> &elements, std::vector<eslocal> &esize, std::vector<eslocal> &etype, const std::vector<eslocal> &indices, CubeFace face) const;

	void pushNodes(std::vector<eslocal> &nodes, const std::vector<eslocal> &indices, CubeFace face) const;
};

}


#endif /* SRC_INPUT_MESHGENERATOR_ELEMENTS_3D_PRISMA15_H_ */
