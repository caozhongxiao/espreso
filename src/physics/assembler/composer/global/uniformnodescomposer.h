
#ifndef SRC_PHYSICS_ASSEMBLER_COMPOSER_GLOBAL_UNIFORMNODESCOMPOSER_H_
#define SRC_PHYSICS_ASSEMBLER_COMPOSER_GLOBAL_UNIFORMNODESCOMPOSER_H_

#include "globalcomposer.h"

namespace espreso {

class UniformNodesComposer: public GlobalComposer {

public:
	UniformNodesComposer(Controller &controler, Provider &provider, int DOFs)
	: GlobalComposer(controler, provider), _DOFs(DOFs) {}

	void initDOFs();
	void buildDirichlet();
	void buildPatterns();
	void buildMVData();

	void assemble(Matrices matrices, const SolverParameters &parameters);
	void setDirichlet(Matrices matrices, const std::vector<double> &subtraction);

	void fillSolution();

protected:
	void synchronize(Matrices matrices);

	void gather(std::vector<double> &data);

	int _DOFs;
};

}


#endif /* SRC_PHYSICS_ASSEMBLER_COMPOSER_GLOBAL_UNIFORMNODESCOMPOSER_H_ */
