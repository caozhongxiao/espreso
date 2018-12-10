
#ifndef SRC_PHYSICS_ASSEMBLER_COMPOSER_GLOBAL_GLOBALCOMPOSER_H_
#define SRC_PHYSICS_ASSEMBLER_COMPOSER_GLOBAL_GLOBALCOMPOSER_H_

#include "../composer.h"

namespace espreso {

class GlobalComposer: public Composer {

public:

	GlobalComposer(Mesh &mesh, DataHolder &instance, Controler &controler)
	: Composer(mesh, instance, controler), _localKOffset(0), _localRHSOffset(0) {}

protected:
	eslocal _localKOffset, _localRHSOffset;
	std::vector<eslocal> _nKSize, _nRHSSize;
	std::vector<eslocal> _tKOffsets, _tRHSOffsets;
	std::vector<eslocal> _KPermutation, _RHSPermutation;
	std::vector<eslocal> _nDistribution;
};

}



#endif /* SRC_PHYSICS_ASSEMBLER_COMPOSER_GLOBAL_GLOBALCOMPOSER_H_ */
