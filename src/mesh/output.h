
#ifndef SRC_NEWOUTPUT_OUTPUT_H_
#define SRC_NEWOUTPUT_OUTPUT_H_

#include <string>

namespace espreso {

class ElementStore;
class BoundaryStore;
class DomainStore;
class RegionStore;

struct NewOutput {

	static void VTKLegacy(const std::string &name, ElementStore *nodes, RegionStore *region);
	static void VTKLegacy(const std::string &name, ElementStore *nodes, DomainStore *domains);
	static void VTKLegacy(const std::string &name, ElementStore *elements, ElementStore *nodes, DomainStore *domains);
	static void VTKLegacy(const std::string &name, BoundaryStore *elements, ElementStore *nodes, bool inner = false);
	static void VTKLegacyDual(const std::string &name, ElementStore *elements, ElementStore *nodes, DomainStore *domains);
};

}




#endif /* SRC_NEWOUTPUT_OUTPUT_H_ */
