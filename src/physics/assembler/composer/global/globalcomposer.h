
#ifndef SRC_PHYSICS_ASSEMBLER_COMPOSER_GLOBAL_GLOBALCOMPOSER_H_
#define SRC_PHYSICS_ASSEMBLER_COMPOSER_GLOBAL_GLOBALCOMPOSER_H_

#include "physics/assembler/composer/composer.h"

namespace espreso {

class Provider;

class GlobalComposer: public Composer {

public:

	GlobalComposer(Controller &controler, Provider &provider)
	: Composer(controler), _provider(provider), _localKOffset(0), _localRHSOffset(0), _MVValuesOffset(0) {}

	void buildMVData();
	void assemble(Matrices matrices, const SolverParameters &parameters);
	void setDirichlet(Matrices matrices, double reduction, const std::vector<double> &subtraction);
	void fillSolution();

	void KplusAlfaM(double alfa);
	void alfaKplusBetaM(double alfa, double beta);
	void enrichRHS(double alfa, NodeData* x);
	void computeReactionForces();
	double residualNormNumerator();
	double residualNormDenominator();

protected:
	virtual void synchronize(Matrices matrices) =0;
	void apply(std::vector<SparseMatrix> &matrices, std::vector<double> &result, std::vector<double> &x);
	virtual void gather(std::vector<double> &data) =0;
	void gather(std::vector<double> &out, std::vector<std::vector<double> > &in)
	{
		out = in[0];
		gather(out);
	}

	Provider &_provider;

	esint _localKOffset, _localRHSOffset;
	std::vector<esint> _nKSize, _nRHSSize;
	std::vector<esint> _tKOffsets, _tRHSOffsets;
	std::vector<esint> _KPermutation, _RHSPermutation;
	std::vector<esint> _nDistribution;

	// data for Mat-Vec product
	esint _MVValuesOffset;
	std::vector<esint> _MVRows, _MVCols;
	std::vector<double> _MVVec;
	std::vector<std::vector<esint> > _MVSend, _MVRecv;
	std::vector<int> _MVNeighbours;

	// data for dirichlet
	std::vector<double> _KDirichletValues;
};

}



#endif /* SRC_PHYSICS_ASSEMBLER_COMPOSER_GLOBAL_GLOBALCOMPOSER_H_ */
