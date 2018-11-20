
#ifndef SRC_PHYSICS_ASSEMBLER_COMPOSER_GLOBAL_GLOBALCOMPOSER_H_
#define SRC_PHYSICS_ASSEMBLER_COMPOSER_GLOBAL_GLOBALCOMPOSER_H_

#include "../composer.h"

namespace espreso {

class GlobalComposer: public Composer {

public:

	GlobalComposer(Mesh &mesh, Step &step, Instance &instance, Controler &controler)
	: Composer(mesh, step, instance, controler) {}
};

}



#endif /* SRC_PHYSICS_ASSEMBLER_COMPOSER_GLOBAL_GLOBALCOMPOSER_H_ */
