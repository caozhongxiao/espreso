
#ifndef SRC_ASSEMBLER_PHYSICS_HEATTRANSFER_H_
#define SRC_ASSEMBLER_PHYSICS_HEATTRANSFER_H_

#include "physics.h"

namespace espreso {

struct HeatTransferConfiguration;
struct ConvectionConfiguration;
struct ResultsSelectionConfiguration;
struct Point;

struct HeatTransfer: public virtual Physics
{
	HeatTransfer(const HeatTransferConfiguration &configuration, const ResultsSelectionConfiguration &propertiesConfiguration);

	virtual MatrixType getMatrixType(const Step &step, size_t domain) const;
	virtual void prepare();
	virtual void preprocessData(const Step &step);
	virtual void analyticRegularization(size_t domain, bool ortogonalCluster);

	virtual ~HeatTransfer() {}

protected:
	void computeInitialTemperature(const Step &step, std::vector<std::vector<double> > &data);

	double computeHTC(const Step &step, const ConvectionConfiguration *convection, const Point &p, double temp) const;

	void convectionMatParameters(
			const ConvectionConfiguration &convection, eslocal eindex, const Point &p, Step step,
			double temp, double T_EXT,
			double &rho, double &dynamic_viscosity, double &dynamic_viscosity_T, double &heat_capacity, double &thermal_conductivity) const;

	const HeatTransferConfiguration &_configuration;
	const ResultsSelectionConfiguration &_propertiesConfiguration;

	NodeData *_temperature;
	ElementData *_gradient, *_flux;
	NodeData *_phaseChange, *_latentHeat;
};

}



#endif /* SRC_ASSEMBLER_PHYSICS_HEATTRANSFER_H_ */
