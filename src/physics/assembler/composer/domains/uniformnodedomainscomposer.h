
#ifndef SRC_PHYSICS_ASSEMBLER_COMPOSER_DOMAINS_UNIFORMNODEDOMAINSCOMPOSER_H_
#define SRC_PHYSICS_ASSEMBLER_COMPOSER_DOMAINS_UNIFORMNODEDOMAINSCOMPOSER_H_

#include "domainscomposer.h"

namespace espreso {

class UniformNodeDomainsComposer: public DomainsComposer {

public:
	UniformNodeDomainsComposer(Mesh &mesh, Step &step, Instance &instance, Controler &controler, int DOFs)
	: DomainsComposer(mesh, step, instance, controler), _DOFs(DOFs) {}

	void initDOFs();
	void buildPatterns();

	void assemble(Matrices matrices);

	void setDirichlet();
	void synchronize();

protected:
	eslocal _DOFs;
};

}


#endif /* SRC_PHYSICS_ASSEMBLER_COMPOSER_DOMAINS_UNIFORMNODEDOMAINSCOMPOSER_H_ */
