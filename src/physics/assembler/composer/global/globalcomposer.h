
#ifndef SRC_PHYSICS_ASSEMBLER_COMPOSER_GLOBAL_GLOBALCOMPOSER_H_
#define SRC_PHYSICS_ASSEMBLER_COMPOSER_GLOBAL_GLOBALCOMPOSER_H_

#include "physics/assembler/composer/composer.h"

namespace espreso {

class Provider;

class GlobalComposer: public Composer {

public:

	GlobalComposer(Controler &controler, Provider &provider)
	: Composer(controler), _provider(provider), _localKOffset(0), _localRHSOffset(0) {}

	NodeData* RHS();

	void KplusAlfaM(double alfa);
	void enrichRHS(double alfa, NodeData* a);
	void RHSMinusR();
	void DirichletMinusRHS();
	double residualNorm();

protected:
	void apply(std::vector<SparseMatrix> &matrices, NodeData *result, NodeData *x);

	Provider &_provider;

	esint _localKOffset, _localRHSOffset;
	std::vector<esint> _nKSize, _nRHSSize;
	std::vector<esint> _tKOffsets, _tRHSOffsets;
	std::vector<esint> _KPermutation, _RHSPermutation;
	std::vector<esint> _nDistribution;
};

}



#endif /* SRC_PHYSICS_ASSEMBLER_COMPOSER_GLOBAL_GLOBALCOMPOSER_H_ */
