
#ifndef SRC_APP_FACTORY_HEATTRANSFERFACTORY_H_
#define SRC_APP_FACTORY_HEATTRANSFERFACTORY_H_

#include "factory.h"

namespace espreso {

struct HeatTransferConfiguration;

class HeatTransferFactory: public FactoryLoader {

public:
	HeatTransferFactory(const HeatTransferConfiguration &configuration, Mesh *mesh);

	size_t loadSteps() const;
	LoadStepSolver* getLoadStepSolver(size_t step, Mesh *mesh, Store *store);

protected:
	const HeatTransferConfiguration &_configuration;
	bool _bem;

};

}



#endif /* SRC_APP_FACTORY_HEATTRANSFERFACTORY_H_ */
