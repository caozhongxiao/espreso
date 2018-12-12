
#ifndef SRC_PHYSICS_ASSEMBLER_COMPOSER_FETI_FETICOMPOSER_H_
#define SRC_PHYSICS_ASSEMBLER_COMPOSER_FETI_FETICOMPOSER_H_

#include "../composer.h"

#include <cstddef>
#include <vector>

namespace espreso {

template <typename TEBoundaries, typename TEData> class serializededata;

class FETIComposer: public Composer {

public:
	FETIComposer(Mesh &mesh, DataHolder &instance, Controler &controler)
	: Composer(mesh, instance, controler) {}

protected:
	std::vector<std::vector<eslocal> > _KPermutation, _RHSPermutation;
	std::vector<size_t> _domainDOFsSize, _domainDirichletSize;
};

}


#endif /* SRC_PHYSICS_ASSEMBLER_COMPOSER_FETI_FETICOMPOSER_H_ */
