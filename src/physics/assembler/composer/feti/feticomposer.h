
#ifndef SRC_PHYSICS_ASSEMBLER_COMPOSER_FETI_FETICOMPOSER_H_
#define SRC_PHYSICS_ASSEMBLER_COMPOSER_FETI_FETICOMPOSER_H_

#include "../composer.h"

#include <cstddef>
#include <vector>

namespace espreso {

class FETIProvider;

template <typename TEBoundaries, typename TEData> class serializededata;
struct FETISolverConfiguration;

class FETIComposer: public Composer {

public:
	FETIComposer(Controler &controler, FETIProvider &provider, FETISolverConfiguration &configuration)
	: Composer(controler), _provider(provider), _configuration(configuration) {}

	NodeData* RHS();

	void KplusAlfaM(double alfa);
	void applyM(NodeData *y, NodeData *x);
	void applyOriginalK(NodeData *y, NodeData *x);
	void enrichRHS(double alfa, NodeData* a);
	void RHSMinusR();
	void DirichletMinusRHS();
	double residualNorm();

protected:
	FETIProvider &_provider;
	FETISolverConfiguration &_configuration;

	std::vector<std::vector<esint> > _KPermutation, _RHSPermutation;
	std::vector<size_t> _domainDOFsSize, _domainDirichletSize;
};

}


#endif /* SRC_PHYSICS_ASSEMBLER_COMPOSER_FETI_FETICOMPOSER_H_ */
