
#ifndef SRC_GLOBALS_RUN_H_
#define SRC_GLOBALS_RUN_H_

#include "../config/ecf/root.h"
#include "../mesh/mesh.h"
#include "../physics/dataholder.h"

namespace espreso {

struct run {

	static ECFRoot ecf;
	static Mesh mesh;
	static DataHolder data;
};

}

#endif /* SRC_GLOBALS_RUN_H_ */
