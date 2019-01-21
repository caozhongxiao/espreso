
#ifndef SRC_PHYSICS_ASSEMBLER_COMPOSER_GLOBAL_UNIFORMNODESCOMPOSER_H_
#define SRC_PHYSICS_ASSEMBLER_COMPOSER_GLOBAL_UNIFORMNODESCOMPOSER_H_

#include "globalcomposer.h"

namespace espreso {

class UniformNodesComposer: public GlobalComposer {

public:
	UniformNodesComposer(Controler &controler, Provider &provider, int DOFs)
	: GlobalComposer(controler, provider), _DOFs(DOFs) {}

	void initDOFs();
	void buildDirichlet();
	void buildPatterns();

	void assemble(Matrices matrices, const SolverParameters &parameters);

	void setDirichlet();
	void synchronize();

	void fillSolution();

protected:
	int _DOFs;
};

}


#endif /* SRC_PHYSICS_ASSEMBLER_COMPOSER_GLOBAL_UNIFORMNODESCOMPOSER_H_ */
