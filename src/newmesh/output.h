
#ifndef SRC_NEWOUTPUT_OUTPUT_H_
#define SRC_NEWOUTPUT_OUTPUT_H_

#include <string>

namespace espreso {

class ElementStore;
class BoundaryStore;

struct NewOutput {

	static void VTKLegacy(const std::string &name, ElementStore *elements, ElementStore *nodes);
	static void VTKLegacy(const std::string &name, BoundaryStore *elements, ElementStore *nodes, bool inner = false);
};

}




#endif /* SRC_NEWOUTPUT_OUTPUT_H_ */
