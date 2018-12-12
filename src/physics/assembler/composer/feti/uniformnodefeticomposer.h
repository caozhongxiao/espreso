
#ifndef SRC_PHYSICS_ASSEMBLER_COMPOSER_FETI_UNIFORMNODEFETICOMPOSER_H_
#define SRC_PHYSICS_ASSEMBLER_COMPOSER_FETI_UNIFORMNODEFETICOMPOSER_H_

#include "../feti/feticomposer.h"

namespace espreso {

class UniformNodeFETIComposer: public FETIComposer {

public:
	UniformNodeFETIComposer(Mesh &mesh, DataHolder &instance, Controler &controler, int DOFs, bool redundantLagrange, bool scaling)
	: FETIComposer(mesh, instance, controler), _DOFs(DOFs), _redundantLagrange(redundantLagrange), _scaling(scaling) {}

	void initDOFs();
	void initDirichlet();
	void buildPatterns();

	void assemble(Matrices matrices);

	void setDirichlet();
	void synchronize();

	void fillSolution();

protected:
	void buildKPattern();
	void buildB1Pattern();
	void buildB0Pattern();
	void updateDuplicity();

	eslocal _DOFs;
	bool _redundantLagrange, _scaling;
};

}


#endif /* SRC_PHYSICS_ASSEMBLER_COMPOSER_FETI_UNIFORMNODEFETICOMPOSER_H_ */
