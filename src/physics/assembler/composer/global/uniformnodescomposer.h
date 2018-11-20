
#ifndef SRC_PHYSICS_ASSEMBLER_COMPOSER_GLOBAL_UNIFORMNODESCOMPOSER_H_
#define SRC_PHYSICS_ASSEMBLER_COMPOSER_GLOBAL_UNIFORMNODESCOMPOSER_H_

#include "globalcomposer.h"

namespace espreso {

class UniformNodesComposer: public GlobalComposer {

public:
	UniformNodesComposer(Mesh &mesh, Step &step, Instance &instance, Controler &controler, int DOFs)
	: GlobalComposer(mesh, step, instance, controler), _DOFs(DOFs) {}

	void initDOFs();
	void buildPatterns();

	void assemble(Matrices matrices);

	void setDirichlet();
	void synchronize();

protected:
	int _DOFs;
};

}


#endif /* SRC_PHYSICS_ASSEMBLER_COMPOSER_GLOBAL_UNIFORMNODESCOMPOSER_H_ */
