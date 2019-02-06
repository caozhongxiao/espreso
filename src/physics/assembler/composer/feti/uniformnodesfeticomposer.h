
#ifndef SRC_PHYSICS_ASSEMBLER_COMPOSER_FETI_UNIFORMNODESFETICOMPOSER_H_
#define SRC_PHYSICS_ASSEMBLER_COMPOSER_FETI_UNIFORMNODESFETICOMPOSER_H_

#include "feticomposer.h"

namespace espreso {

class UniformNodesFETIComposer: public FETIComposer {

public:
	UniformNodesFETIComposer(Controller &controler, FETIProvider &provider, FETISolverConfiguration &configuration, int DOFs)
	: FETIComposer(controler, provider, configuration), _DOFs(DOFs) {}

	void initDOFs();
	void buildDirichlet();
	void buildPatterns();
	void buildMVData() {} // no data are needed

	void assemble(Matrices matrices, const SolverParameters &parameters);
	void setDirichlet(Matrices matrices, double reduction, const std::vector<double> &subtraction);

	void fillSolution();

protected:
	void buildKFEMPattern(esint domain);
	void buildKBEMPattern(esint domain);

	void buildB1Pattern();
	void updateDuplicity();

	void divide(std::vector<double> &in, std::vector<std::vector<double> > &out)
	{
		_divide(in, out, true);
	}
	void duply(std::vector<double> &in, std::vector<std::vector<double> > &out)
	{
		_divide(in, out, false);
	}
	void avgGather(std::vector<double> &out, std::vector<std::vector<double> > &in)
	{
		_gather(out, in, true);
	}
	void gather(std::vector<double> &out, std::vector<std::vector<double> > &in)
	{
		_gather(out, in, false);
	}

	void _divide(std::vector<double> &in, std::vector<std::vector<double> > &out, bool split);
	void _gather(std::vector<double> &out, std::vector<std::vector<double> > &in, bool avg);

	esint _DOFs;
};

}


#endif /* SRC_PHYSICS_ASSEMBLER_COMPOSER_FETI_UNIFORMNODESFETICOMPOSER_H_ */
