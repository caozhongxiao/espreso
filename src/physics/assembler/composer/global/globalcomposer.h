
#ifndef SRC_PHYSICS_ASSEMBLER_COMPOSER_GLOBAL_GLOBALCOMPOSER_H_
#define SRC_PHYSICS_ASSEMBLER_COMPOSER_GLOBAL_GLOBALCOMPOSER_H_

#include "physics/assembler/composer/composer.h"

namespace espreso {

class Provider;

class GlobalComposer: public Composer {

public:

	GlobalComposer(Controller &controler, Provider &provider)
	: Composer(controler), _provider(provider), _localKOffset(0), _localRHSOffset(0) {}

	void KplusAlfaM(double alfa);
	void computeReactionForces();
	double residualNormNumerator();
	double residualNormDenominator();

protected:
	void apply(std::vector<SparseMatrix> &matrices, std::vector<double> &result, std::vector<double> &x);
	virtual void gather(std::vector<double> &data) =0;

	Provider &_provider;

	esint _localKOffset, _localRHSOffset;
	std::vector<esint> _nKSize, _nRHSSize;
	std::vector<esint> _tKOffsets, _tRHSOffsets;
	std::vector<esint> _KPermutation, _RHSPermutation;
	std::vector<esint> _nDistribution;

	// data for Mat-Vec product
	std::vector<esint> _MVCols;
	std::vector<double> _MVVec;
	std::vector<std::vector<esint> > _MVSend, _MVRecv;
};

}



#endif /* SRC_PHYSICS_ASSEMBLER_COMPOSER_GLOBAL_GLOBALCOMPOSER_H_ */
