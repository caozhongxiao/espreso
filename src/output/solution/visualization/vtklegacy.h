
#ifndef SRC_OUTPUT_SOLUTION_VISUALIZATION_VTKLEGACY_H_
#define SRC_OUTPUT_SOLUTION_VISUALIZATION_VTKLEGACY_H_

#include "visualization.h"

#include <string>

namespace espreso {

struct OutputConfiguration;
class ElementStore;
class NodeStore;

struct VTKLegacy: public Visualization {

	static void mesh(const std::string &name, const OutputConfiguration &configuration, NodeStore *nodes, ElementStore *elements);
	static void solution(const std::string &name, const OutputConfiguration &configuration, NodeStore *nodes, ElementStore *elements);
	static void nodesIntervals(const std::string &name, NodeStore *nodes);

//	static void VTKLegacy(const std::string &name, ElementStore *nodes, RegionStore *region);
//	static void VTKLegacy(const std::string &name, ElementStore *nodes, DomainStore *domains);
//	static void VTKLegacy(const std::string &name, ElementStore *elements, ElementStore *nodes, DomainStore *domains);
//	static void VTKLegacy(const std::string &name, BoundaryStore *elements, ElementStore *nodes, bool inner = false);
//	static void VTKLegacyDual(const std::string &name, ElementStore *elements, ElementStore *nodes, DomainStore *domains);
};

}


#endif /* SRC_OUTPUT_SOLUTION_VISUALIZATION_VTKLEGACY_H_ */
