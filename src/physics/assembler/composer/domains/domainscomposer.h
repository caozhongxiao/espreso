
#ifndef SRC_PHYSICS_ASSEMBLER_COMPOSER_DOMAINS_DOMAINSCOMPOSER_H_
#define SRC_PHYSICS_ASSEMBLER_COMPOSER_DOMAINS_DOMAINSCOMPOSER_H_

#include "../composer.h"

#include <cstddef>
#include <vector>

namespace espreso {

template <typename TEBoundaries, typename TEData> class serializededata;

class DomainsComposer: public Composer {

public:
	DomainsComposer(Mesh &mesh, DataHolder &instance, Controler &controler)
	: Composer(mesh, instance, controler) {}

protected:
	std::vector<std::vector<eslocal> > _KPermutation, _RHSPermutation;
	std::vector<size_t> _domainDirichletSize;
};

}


#endif /* SRC_PHYSICS_ASSEMBLER_COMPOSER_DOMAINS_DOMAINSCOMPOSER_H_ */
