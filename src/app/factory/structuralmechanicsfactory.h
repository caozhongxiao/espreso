
#ifndef SRC_APP_FACTORY_STRUCTURALMECHANICSFACTORY_H_
#define SRC_APP_FACTORY_STRUCTURALMECHANICSFACTORY_H_

#include "factory.h"

namespace espreso {

struct StructuralMechanicsConfiguration;
struct ResultsSelectionConfiguration;

class StructuralMechanicsFactory: public FactoryLoader {

public:
	StructuralMechanicsFactory(StructuralMechanicsConfiguration &configuration, ResultsSelectionConfiguration &propertiesConfiguration, Mesh *mesh);

	size_t loadSteps() const;
	LoadStepSolver* getLoadStepSolver(size_t step, Mesh *mesh);

protected:
	StructuralMechanicsConfiguration &_configuration;
	ResultsSelectionConfiguration &_propertiesConfiguration;
	bool _bem;

};

}


#endif /* SRC_APP_FACTORY_STRUCTURALMECHANICSFACTORY_H_ */
