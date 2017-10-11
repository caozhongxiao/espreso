
#ifndef SRC_NEWOUTPUT_OUTPUT_H_
#define SRC_NEWOUTPUT_OUTPUT_H_

namespace espreso {

class ElementStore;

class NewOutput {

	static void VTKLegacy(ElementStore *elements, ElementStore *nodes);
};

}




#endif /* SRC_NEWOUTPUT_OUTPUT_H_ */
