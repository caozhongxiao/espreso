
#ifndef SRC_ASSEMBLER_PHYSICS_HEATTRANSFER_H_
#define SRC_ASSEMBLER_PHYSICS_HEATTRANSFER_H_

#include "../assembler/physics.h"

namespace espreso {

struct HeatTransferConfiguration;
struct ConvectionConfiguration;
struct ResultsSelectionConfiguration;
struct Point;
struct BoundaryRegionStore;

struct HeatTransfer: public virtual Physics
{
	HeatTransfer(const HeatTransferConfiguration &configuration, const ResultsSelectionConfiguration &propertiesConfiguration);

	virtual void initLocalDOFs(std::vector<eslocal> &offsets) { initLocalNodeUniformDOFs(offsets, 1); }
	virtual void initGlobalDOFs(eslocal &offset) { initGlobalNodeUniformDOFs(offset, 1); }

	virtual void buildLocalCSRPattern() { buildLocalNodeUniformCSRPattern(1); }
	virtual void buildGlobalCSRPattern() { buildGlobalNodeUniformCSRPattern(1); }

	virtual MatrixType getMatrixType(size_t domain) const;
	virtual void prepare();
	virtual void preprocessData();
	virtual void setDirichlet();
	virtual void analyticRegularization(size_t domain, bool ortogonalCluster);

	double sumSquares(const std::vector<std::vector<double> > &data, SumRestriction restriction) const;

	virtual ~HeatTransfer() {}

protected:
	void computeInitialTemperature(std::vector<std::vector<double> > &data);

	double computeHTC(const ConvectionConfiguration *convection, const Point &p, double temp) const;

	void convectionMatParameters(
			const ConvectionConfiguration *convection, const Point &p, double temp, double T_EXT,
			double &rho, double &dynamic_viscosity, double &dynamic_viscosity_T, double &heat_capacity, double &thermal_conductivity) const;

	const HeatTransferConfiguration &_configuration;
	const ResultsSelectionConfiguration &_propertiesConfiguration;

	NodeData *_temperature;
	ElementData *_gradient, *_flux;
	NodeData *_phaseChange, *_latentHeat;
	ElementData *_translationMotion;
};

}



#endif /* SRC_ASSEMBLER_PHYSICS_HEATTRANSFER_H_ */
