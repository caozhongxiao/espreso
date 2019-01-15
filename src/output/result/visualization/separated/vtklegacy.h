
#ifndef SRC_OUTPUT_RESULT_VISUALIZATION_SEPARATED_VTKLEGACY_H_
#define SRC_OUTPUT_RESULT_VISUALIZATION_SEPARATED_VTKLEGACY_H_

#include <string>
#include <vector>

#include "separatedvisualization.h"

namespace espreso {

struct DataHolder;
struct SpaceFillingCurve;

struct VTKLegacy: public SeparatedVisualization {

	static double clusterShrinkRatio, domainShrinkRatio;

protected:
	VTKLegacy(const Mesh &mesh);

	void mesh(const std::string &name);
	void solution(const std::string &name);
	void nodesIntervals(const std::string &name);
	void externalIntervals(const std::string &name);
	void sharedInterface(const std::string &name);
	void surface(const std::string &name);
	void domainSurface(const std::string &name);
	void corners(const std::string &name);
	void sFixPoints(const std::string &name);
	void iFixPoints(const std::string &name);
	void contact(const std::string &name);
	void closeElements(const std::string &name);
	void neighbors(const std::string &name);

	void points(const std::string &name, const std::vector<esint> &points, const std::vector<esint> &distribution);
};

struct VTKLegacyDebugInfo: public VTKLegacy {

	VTKLegacyDebugInfo(const Mesh &mesh);

	static void dirichlet(const Mesh &mesh, const DataHolder &instance);
	static void gluing(const Mesh &mesh, const DataHolder &instance);

	static void spaceFillingCurve(const SpaceFillingCurve &sfc, const std::vector<uint> &bucketsBorders);

	void updateMesh()
	{
		mesh(_path + "mesh");
		nodesIntervals(_path + "nodeintervals");
		externalIntervals(_path + "externalIntervals");
		sharedInterface(_path + "sharedinterfaces");
		surface(_path + "surface");
		domainSurface(_path + "domainSurface");
		corners(_path + "corners");
		sFixPoints(_path + "sfixpoints");
		iFixPoints(_path + "ifixpoints");
		contact(_path + "contact");
		closeElements(_path + "closeElements");
		neighbors(_path + "eneighbors");
	}
	void updateSolution()
	{
		solution(_path + "solution");
	}

protected:
	std::string _path;
};

}


#endif /* SRC_OUTPUT_RESULT_VISUALIZATION_SEPARATED_VTKLEGACY_H_ */
