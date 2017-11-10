
#ifndef SRC_APP_FACTORY_STRUCTURALMECHANICSFACTORY_H_
#define SRC_APP_FACTORY_STRUCTURALMECHANICSFACTORY_H_

#include "factory.h"

namespace espreso {

struct StructuralMechanicsConfiguration;
struct ResultsSelectionConfiguration;

class StructuralMechanicsFactory: public FactoryLoader {

public:
	StructuralMechanicsFactory(const StructuralMechanicsConfiguration &configuration, const ResultsSelectionConfiguration &propertiesConfiguration, OldMesh *mesh);

	size_t loadSteps() const;
	LoadStepSolver* getLoadStepSolver(size_t step, OldMesh *mesh, ResultStoreList *store);

protected:
	const StructuralMechanicsConfiguration &_configuration;
	const ResultsSelectionConfiguration &_propertiesConfiguration;
	bool _bem;

};

}


#endif /* SRC_APP_FACTORY_STRUCTURALMECHANICSFACTORY_H_ */
