
#ifndef SRC_APP_FACTORY_STRUCTURALMECHANICSFACTORY_H_
#define SRC_APP_FACTORY_STRUCTURALMECHANICSFACTORY_H_

#include "factory.h"

namespace espreso {

struct StructuralMechanicsConfiguration;

class StructuralMechanicsFactory: public FactoryLoader {

public:
	StructuralMechanicsFactory(const StructuralMechanicsConfiguration &configuration, Mesh *mesh);

	size_t loadSteps() const;
	LoadStepSolver* getLoadStepSolver(size_t step, Mesh *mesh, Store *store);

protected:
	const StructuralMechanicsConfiguration &_configuration;
	bool _bem;

};

}


#endif /* SRC_APP_FACTORY_STRUCTURALMECHANICSFACTORY_H_ */
