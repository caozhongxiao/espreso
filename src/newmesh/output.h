
#ifndef SRC_NEWOUTPUT_OUTPUT_H_
#define SRC_NEWOUTPUT_OUTPUT_H_

#include <string>

namespace espreso {

class ElementStore;

struct NewOutput {

	static void VTKLegacy(const std::string &name, ElementStore *elements, ElementStore *nodes);
};

}




#endif /* SRC_NEWOUTPUT_OUTPUT_H_ */
