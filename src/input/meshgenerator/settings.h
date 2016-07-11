
#ifndef INPUT_MESHGENERATOR_SETTINGS_H_
#define INPUT_MESHGENERATOR_SETTINGS_H_

#include "esbasis.h"
#include "esmesh.h"

namespace espreso {
namespace input {

enum class GeneratorShape {
	CUBE,
	SPHERE,
	PLANE
};

enum class ElementType {
	HEXA8,
	HEXA20,
	TETRA4,
	TETRA10,
	PRISMA6,
	PRISMA15,
	PYRAMID5,
	PYRAMID13,

	SQUARE4,
	SQUARE8,
	TRIANGLE3,
	TRIANGLE6
};

struct Settings {

	Settings(const Configuration &configuration, size_t index, size_t size);
	Settings(size_t index, size_t size);

	std::vector<Parameter> parameters;

	size_t index;
	size_t size;

	ElementType eType;
	GeneratorShape shape;
	config::assembler::PHYSICSalternative physics;

	std::vector<Material> materials;

	bool useMetis;

protected:
	void defaultSettings();
};

inline std::ostream& operator<<(std::ostream& os, const Settings &s)
{
       os << "clusters: " << s.size << "\n";
       os << "generated shape: " << static_cast<int>(s.shape) << "\n";
       os << "type of the element: " << static_cast<int>(s.eType) << "\n";
       os << "physics: " << static_cast<int>(s.physics) << "\n";

       os << "material 1:\n";
       os << "\tdensity: " << s.materials[0].density << "\n";
       os << "\tyoung's modulus: " << s.materials[0].youngModulus << "\n";
       os << "\tpoisson's ratio: " << s.materials[0].poissonRatio << "\n";

       os << "material 2:\n";
       os << "\tdensity: " << s.materials[1].density << "\n";
       os << "\tyoung's modulus: " << s.materials[1].youngModulus << "\n";
       os << "\tpoisson's ratio: " << s.materials[1].poissonRatio << "\n";

       os << "partition: " << (s.useMetis ? "by METIS" : "regular") << "\n";
       return os;
}

}
}



#endif /* INPUT_MESHGENERATOR_SETTINGS_H_ */

