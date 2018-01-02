
#ifndef SRC_OUTPUT_RESULT_VISUALIZATION_VTKLEGACY_H_
#define SRC_OUTPUT_RESULT_VISUALIZATION_VTKLEGACY_H_

#include <string>

#include "visualization.h"

namespace espreso {

struct VTKLegacy: public Visualization {

protected:
	VTKLegacy(const Mesh &mesh, double clusterShrinkRatio, double domainShrinkRatio);

	void mesh(const std::string &name);
	void solution(const std::string &name);
	void nodesIntervals(const std::string &name);
	void externalIntervals(const std::string &name);
	void sharedInterface(const std::string &name);
	void corners(const std::string &name);

	double _clusterShrinkRatio, _domainShrinkRatio;
};

struct VTKLegacyDebugInfo: public VTKLegacy {

	VTKLegacyDebugInfo(const Mesh &mesh, double clusterShrinkRatio, double domainShrinkRatio);

	void updateMesh()
	{
		mesh(_path + "mesh");
		nodesIntervals(_path + "nodeintervals");
		externalIntervals(_path + "externalIntervals");
		sharedInterface(_path + "sharedinterfaces");
		corners(_path + "corners");
	}
	void updateSolution(const Step &step)
	{
		solution(_path + "solution");
	}

protected:
	std::string _path;
};

}


#endif /* SRC_OUTPUT_RESULT_VISUALIZATION_VTKLEGACY_H_ */
