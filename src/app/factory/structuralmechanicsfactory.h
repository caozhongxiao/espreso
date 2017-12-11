
#ifndef SRC_APP_FACTORY_STRUCTURALMECHANICSFACTORY_H_
#define SRC_APP_FACTORY_STRUCTURALMECHANICSFACTORY_H_

#include "factory.h"

namespace espreso {

struct StructuralMechanicsConfiguration;
struct ResultsSelectionConfiguration;

class StructuralMechanicsFactory: public FactoryLoader {

public:
	StructuralMechanicsFactory(const StructuralMechanicsConfiguration &configuration, const ResultsSelectionConfiguration &propertiesConfiguration, Mesh *mesh);

	size_t loadSteps() const;
	LoadStepSolver* getLoadStepSolver(size_t step, Mesh *mesh, ResultStore *store);

protected:
	const StructuralMechanicsConfiguration &_configuration;
	const ResultsSelectionConfiguration &_propertiesConfiguration;
	bool _bem;

};

}


#endif /* SRC_APP_FACTORY_STRUCTURALMECHANICSFACTORY_H_ */
