
#ifndef SRC_INPUT_GENERATOR_GENERATOR_H_
#define SRC_INPUT_GENERATOR_GENERATOR_H_

#include <cstddef>

namespace espreso {

struct InputGeneratorConfiguration;
class OldMesh;

namespace input {

struct Generator {

	static void generate(const InputGeneratorConfiguration &configuration, OldMesh &mesh, size_t index, size_t size);

	static double precision;
};
}
}



#endif /* SRC_INPUT_GENERATOR_GENERATOR_H_ */
