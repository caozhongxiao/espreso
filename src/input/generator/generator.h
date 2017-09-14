
#ifndef SRC_INPUT_GENERATOR_GENERATOR_H_
#define SRC_INPUT_GENERATOR_GENERATOR_H_

#include <cstddef>

namespace espreso {

struct InputGeneratorConfiguration;
class Mesh;

namespace input {

struct Generator {

	static void generate(const InputGeneratorConfiguration &configuration, Mesh &mesh, size_t index, size_t size);

	static double precision;
};
}
}



#endif /* SRC_INPUT_GENERATOR_GENERATOR_H_ */
