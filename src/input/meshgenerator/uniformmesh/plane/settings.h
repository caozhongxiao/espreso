
#ifndef INPUT_MESHGENERATOR_UNIFORMMESH_PLANE_SETTINGS_H_
#define INPUT_MESHGENERATOR_UNIFORMMESH_PLANE_SETTINGS_H_

#include <map>

#include "../cube/settings.h"

namespace espreso {
namespace input {

struct PlaneSettings: public CubeSettings {

	PlaneSettings(const Configuration &configuration, size_t index, size_t size);
	PlaneSettings(size_t index, size_t size);

protected:
	void defaultPlaneSettings();
};

inline std::ostream& operator<<(std::ostream& os, const PlaneSettings &s)
{
	os << UniformSettings(s);
	os << "clusters: " << s.clusters[0] << " : " << s.clusters[1] << "\n";
	os << "cube length: " << s.problemLength[0] << " : " << s.problemLength[1] << "\n";

	std::vector<std::string> properties = { "DIRICHLET", "FORCES" };
	std::vector<std::string> cube_faces = { "FRONT", "REAR", "LEFT", "RIGHT", "TOP", "BOTTOM" };
	std::vector<std::string> axis = { "X", "Y", "Z" };

	for (size_t f = 2; f < cube_faces.size(); f++) {
		for (size_t p = DIRICHLET_X; p <= FORCES_Z; p++) {
			std::string name = properties[p / 3] + "_" + cube_faces[f] + "_" + axis[p % 3];
			if (s.boundaryCondition[f * 6 + p] != std::numeric_limits<double>::infinity()) {
				os << name << ": " << s.boundaryCondition[f * 6 + p] << "\n";
			}
		}
	}

	return os;
}

}
}

#endif /* INPUT_MESHGENERATOR_UNIFORMMESH_PLANE_SETTINGS_H_ */
