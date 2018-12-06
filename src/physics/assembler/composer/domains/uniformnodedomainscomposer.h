
#ifndef SRC_PHYSICS_ASSEMBLER_COMPOSER_DOMAINS_UNIFORMNODEDOMAINSCOMPOSER_H_
#define SRC_PHYSICS_ASSEMBLER_COMPOSER_DOMAINS_UNIFORMNODEDOMAINSCOMPOSER_H_

#include "domainscomposer.h"

namespace espreso {

class UniformNodeDomainsComposer: public DomainsComposer {

public:
	UniformNodeDomainsComposer(Mesh &mesh, Step &step, Instance &instance, Controler &controler, int DOFs, bool redundantLagrange, bool scaling)
	: DomainsComposer(mesh, step, instance, controler), _DOFs(DOFs), _redundantLagrange(redundantLagrange), _scaling(scaling) {}

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


#endif /* SRC_PHYSICS_ASSEMBLER_COMPOSER_DOMAINS_UNIFORMNODEDOMAINSCOMPOSER_H_ */
