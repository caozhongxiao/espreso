
#ifndef SRC_PHYSICS_ASSEMBLER_COMPOSER_GLOBAL_GLOBALCOMPOSER_H_
#define SRC_PHYSICS_ASSEMBLER_COMPOSER_GLOBAL_GLOBALCOMPOSER_H_

#include "../composer.h"

namespace espreso {

class Provider;

class GlobalComposer: public Composer {

public:

	GlobalComposer(Controler &controler, Provider &provider)
	: Composer(controler), _provider(provider), _localKOffset(0), _localRHSOffset(0) {}

	NodeData* RHS();

	void KplusAlfaM(double alfa);
	void applyM(NodeData *y, NodeData *x);
	void applyOriginalK(NodeData *y, NodeData *x);
	void enrichRHS(double alfa, NodeData* a);
	void RHSMinusR();
	void DirichletMinusRHS();
	double residualNorm();

protected:
	Provider &_provider;

	esint _localKOffset, _localRHSOffset;
	std::vector<esint> _nKSize, _nRHSSize;
	std::vector<esint> _tKOffsets, _tRHSOffsets;
	std::vector<esint> _KPermutation, _RHSPermutation;
	std::vector<esint> _nDistribution;
};

}



#endif /* SRC_PHYSICS_ASSEMBLER_COMPOSER_GLOBAL_GLOBALCOMPOSER_H_ */
