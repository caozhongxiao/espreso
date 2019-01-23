
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
	void buildMVData();

	void assemble(Matrices matrices, const SolverParameters &parameters);
	void setDirichlet();

	void fillSolution();

protected:
	void synchronize(Matrices matrices);

	void gather(NodeData *data);

	int _DOFs;
};

}


#endif /* SRC_PHYSICS_ASSEMBLER_COMPOSER_GLOBAL_UNIFORMNODESCOMPOSER_H_ */
