
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
	FETIComposer(Controller &controler, FETIProvider &provider, FETISolverConfiguration &configuration);

	void KplusAlfaM(double alfa);
	void enrichRHS(double alfa, NodeData* x);
	void computeReactionForces() {} // returned by FETI solver
	double residualNormNumerator();
	double residualNormDenominator();

protected:
	bool isBEMDomain(esint domain);

	void apply(std::vector<SparseMatrix> &matrices, std::vector<double> &result, std::vector<double> &x);
	virtual void divide(std::vector<double> &in, std::vector<std::vector<double> > &out) =0;
	virtual void duply(std::vector<double> &in, std::vector<std::vector<double> > &out) =0;
	virtual void avgGather(std::vector<double> &out, std::vector<std::vector<double> > &in) =0;
	virtual void gather(std::vector<double> &out, std::vector<std::vector<double> > &in) =0;

	FETIProvider &_provider;
	FETISolverConfiguration &_configuration;

	std::vector<std::vector<esint> > _KPermutation, _RHSPermutation;
	std::vector<size_t> _domainDOFsSize, _domainDirichletSize;
	std::vector<int> _BEMDomain;
};

}


#endif /* SRC_PHYSICS_ASSEMBLER_COMPOSER_FETI_FETICOMPOSER_H_ */
