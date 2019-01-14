
#ifndef SRC_PHYSICS_ASSEMBLER_COMPOSER_FETI_UNIFORMNODESFETICOMPOSER_H_
#define SRC_PHYSICS_ASSEMBLER_COMPOSER_FETI_UNIFORMNODESFETICOMPOSER_H_

#include "feticomposer.h"

namespace espreso {

class UniformNodesFETIComposer: public FETIComposer {

public:
	UniformNodesFETIComposer(Controler &controler, FETIProvider &provider, FETISolverConfiguration &configuration, int DOFs)
	: FETIComposer(controler, provider, configuration), _DOFs(DOFs) {}

	void initDOFs();
	void initDirichlet();
	void buildPatterns();

	void assemble(Matrices matrices, const SolverParameters &parameters);

	void setDirichlet();
	void synchronize();

	void fillSolution();

protected:
	void buildKPattern();
	void buildB1Pattern();
	void updateDuplicity();

	esint _DOFs;
};

}


#endif /* SRC_PHYSICS_ASSEMBLER_COMPOSER_FETI_UNIFORMNODESFETICOMPOSER_H_ */
