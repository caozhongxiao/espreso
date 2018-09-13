
#ifndef SRC_INPUT_MESHGENERATOR_ELEMENTS_3D_PRISMA6_H_
#define SRC_INPUT_MESHGENERATOR_ELEMENTS_3D_PRISMA6_H_

#include "linearvolume.h"

namespace espreso {

struct Prisma6Generator: public LinearVolumeGenerator {

	Prisma6Generator();

	void pushElements(std::vector<eslocal> &elements, const std::vector<eslocal> &indices) const;
	void pushFace(std::vector<eslocal> &elements, std::vector<eslocal> &esize, std::vector<int> &etype, const std::vector<eslocal> &indices, CubeFace face) const;

	void pushNodes(std::vector<eslocal> &nodes, const std::vector<eslocal> &indices, CubeFace face) const;
};

}

#endif /* SRC_INPUT_MESHGENERATOR_ELEMENTS_3D_PRISMA6_H_ */
