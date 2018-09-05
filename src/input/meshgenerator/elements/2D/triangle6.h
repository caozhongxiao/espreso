
#ifndef INPUT_MESHGENERATOR_ELEMENTS_2D_TRIANGLE6_H_
#define INPUT_MESHGENERATOR_ELEMENTS_2D_TRIANGLE6_H_

#include "quadraticplane.h"

namespace espreso {

struct Triangle6Generator: public QuadraticPlaneGenerator {

	Triangle6Generator();

	void pushElements(std::vector<eslocal> &elements, const std::vector<eslocal> &indices) const;
};

}


#endif /* INPUT_MESHGENERATOR_ELEMENTS_2D_TRIANGLE6_H_ */
