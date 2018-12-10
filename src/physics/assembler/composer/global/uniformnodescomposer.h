
#ifndef SRC_PHYSICS_ASSEMBLER_COMPOSER_GLOBAL_UNIFORMNODESCOMPOSER_H_
#define SRC_PHYSICS_ASSEMBLER_COMPOSER_GLOBAL_UNIFORMNODESCOMPOSER_H_

#include "globalcomposer.h"

namespace espreso {

class UniformNodesComposer: public GlobalComposer {

public:
	UniformNodesComposer(Mesh &mesh, DataHolder &instance, Controler &controler, int DOFs)
	: GlobalComposer(mesh, instance, controler), _DOFs(DOFs) {}

	void initDOFs();
	void initDirichlet();
	void buildPatterns();

	void assemble(Matrices matrices);

	void setDirichlet();
	void synchronize();

	void fillSolution();

protected:
	int _DOFs;
};

}


#endif /* SRC_PHYSICS_ASSEMBLER_COMPOSER_GLOBAL_UNIFORMNODESCOMPOSER_H_ */
