
#ifndef INPUT_MESHGENERATOR_ELEMENTS_2D_SQUARE4_H_
#define INPUT_MESHGENERATOR_ELEMENTS_2D_SQUARE4_H_

#include "linearplane.h"

namespace espreso {

struct Square4Generator: public LinearPlaneGenerator {

	Square4Generator();

	void pushElements(std::vector<eslocal> &elements, const std::vector<eslocal> &indices) const;
};

}


#endif /* INPUT_MESHGENERATOR_ELEMENTS_2D_SQUARE4_H_ */
