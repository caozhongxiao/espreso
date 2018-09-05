
#ifndef SRC_INPUT_MESHGENERATOR_MESHGENERATOR_H_
#define SRC_INPUT_MESHGENERATOR_MESHGENERATOR_H_

#include <cstddef>

namespace espreso {

class InputGeneratorConfiguration;
class Mesh;

class MeshGenerator {

public:
	static double precision;

	static void generate(const InputGeneratorConfiguration &configuration, Mesh &mesh);

protected:
	MeshGenerator(const InputGeneratorConfiguration &configuration, Mesh &mesh);
};

}



#endif /* SRC_INPUT_MESHGENERATOR_MESHGENERATOR_H_ */
