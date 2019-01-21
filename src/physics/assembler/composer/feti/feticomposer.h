
#ifndef SRC_PHYSICS_ASSEMBLER_COMPOSER_FETI_FETICOMPOSER_H_
#define SRC_PHYSICS_ASSEMBLER_COMPOSER_FETI_FETICOMPOSER_H_

#include "physics/assembler/composer/composer.h"

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
	void enrichRHS(double alfa, NodeData* x);
	void RHSMinusR();
	void DirichletMinusRHS();
	double residualNorm();

protected:
	void apply(std::vector<SparseMatrix> &matrices, NodeData *result, NodeData *x);
	virtual void divide(NodeData *in, std::vector<std::vector<double> > &out) =0;
	virtual void gather(NodeData *out, std::vector<std::vector<double> > &in) =0;

	FETIProvider &_provider;
	FETISolverConfiguration &_configuration;

	std::vector<std::vector<esint> > _KPermutation, _RHSPermutation;
	std::vector<size_t> _domainDOFsSize, _domainDirichletSize;
};

}


#endif /* SRC_PHYSICS_ASSEMBLER_COMPOSER_FETI_FETICOMPOSER_H_ */
