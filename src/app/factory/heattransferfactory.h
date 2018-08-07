
#ifndef SRC_APP_FACTORY_HEATTRANSFERFACTORY_H_
#define SRC_APP_FACTORY_HEATTRANSFERFACTORY_H_

#include "factory.h"

namespace espreso {

struct HeatTransferConfiguration;
struct ResultsSelectionConfiguration;

class HeatTransferFactory: public FactoryLoader {

public:
	HeatTransferFactory(Step *step, const HeatTransferConfiguration &configuration, const ResultsSelectionConfiguration &propertiesConfiguration, Mesh *mesh);

	size_t loadSteps() const;
	LoadStepSolver* getLoadStepSolver(size_t step, Mesh *mesh, ResultStore *store);

protected:
	Step *_step;
	const HeatTransferConfiguration &_configuration;
	const ResultsSelectionConfiguration &_propertiesConfiguration;
	bool _bem;

};

}



#endif /* SRC_APP_FACTORY_HEATTRANSFERFACTORY_H_ */
