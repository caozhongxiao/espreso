
#ifndef INPUT_MESHGENERATOR_ELEMENTS_2D_TRIANGLE3_H_
#define INPUT_MESHGENERATOR_ELEMENTS_2D_TRIANGLE3_H_

#include "linearplane.h"

namespace espreso {

struct Triangle3Generator: public LinearPlaneGenerator {

	Triangle3Generator();

	void pushElements(std::vector<eslocal> &elements, const std::vector<eslocal> &indices) const;
};

}


#endif /* INPUT_MESHGENERATOR_ELEMENTS_2D_TRIANGLE3_H_ */
