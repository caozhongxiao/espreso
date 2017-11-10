
#ifndef SRC_APP_FACTORY_HEATTRANSFERFACTORY_H_
#define SRC_APP_FACTORY_HEATTRANSFERFACTORY_H_

#include "factory.h"

namespace espreso {

struct HeatTransferConfiguration;
struct ResultsSelectionConfiguration;

class HeatTransferFactory: public FactoryLoader {

public:
	HeatTransferFactory(const HeatTransferConfiguration &configuration, const ResultsSelectionConfiguration &propertiesConfiguration, OldMesh *mesh);

	size_t loadSteps() const;
	LoadStepSolver* getLoadStepSolver(size_t step, OldMesh *mesh, ResultStoreList *store);

protected:
	const HeatTransferConfiguration &_configuration;
	const ResultsSelectionConfiguration &_propertiesConfiguration;
	bool _bem;

};

}



#endif /* SRC_APP_FACTORY_HEATTRANSFERFACTORY_H_ */
