
#ifndef INPUT_MESHGENERATOR_ELEMENTS_2D_SQUARE8_H_
#define INPUT_MESHGENERATOR_ELEMENTS_2D_SQUARE8_H_

#include "quadraticplane.h"

namespace espreso {

struct Square8Generator: public QuadraticPlaneGenerator {

	Square8Generator();

	void pushElements(std::vector<eslocal> &elements, const std::vector<eslocal> &indices) const;
};

}


#endif /* INPUT_MESHGENERATOR_ELEMENTS_2D_SQUARE8_H_ */
