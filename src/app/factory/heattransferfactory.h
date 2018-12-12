
#ifndef SRC_APP_FACTORY_HEATTRANSFERFACTORY_H_
#define SRC_APP_FACTORY_HEATTRANSFERFACTORY_H_

#include "factory.h"

namespace espreso {

struct HeatTransferConfiguration;
struct ResultsSelectionConfiguration;

class HeatTransferFactory: public FactoryLoader {

public:
	HeatTransferFactory(HeatTransferConfiguration &configuration, ResultsSelectionConfiguration &propertiesConfiguration, Mesh *mesh);

	size_t loadSteps() const;
	LoadStepSolver* getLoadStepSolver(size_t step, Mesh *mesh);

protected:
	HeatTransferConfiguration &_configuration;
	ResultsSelectionConfiguration &_propertiesConfiguration;
	bool _bem;

};

}



#endif /* SRC_APP_FACTORY_HEATTRANSFERFACTORY_H_ */
