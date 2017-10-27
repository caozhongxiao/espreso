
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

	virtual std::vector<size_t> solutionsIndicesToStore() const;

	virtual MatrixType getMatrixType(const Step &step, size_t domain) const;
	virtual bool isMatrixTimeDependent(const Step &step) const;
	virtual bool isMatrixTemperatureDependent(const Step &step) const;
	virtual void prepare();
	virtual void preprocessData(const Step &step);
	virtual void analyticRegularization(size_t domain, bool ortogonalCluster);

	const std::vector<Property>& pointDOFs() const
	{
		static std::vector<Property> pointDOFs = { Property::TEMPERATURE };
		return pointDOFs;
	}
	const std::vector<Property>& midPointDOFs() const
	{
		static std::vector<Property> midPointDOFs = { Property::TEMPERATURE };
		return midPointDOFs;
	}
	const std::vector<Property>& edgeDOFs() const
	{
		static std::vector<Property> edgeDOFs = { };
		return edgeDOFs;
	}
	const std::vector<Property>& faceDOFs() const
	{
		static std::vector<Property> faceDOFs = { };
		return faceDOFs;
	}
	const std::vector<Property>& elementDOFs() const
	{
		static std::vector<Property> elementDOFs = { };
		return elementDOFs;
	}

	virtual ~HeatTransfer() {}

protected:
	enum SolutionIndex: size_t {
		TEMPERATURE,
		GRADIENT,
		FLUX,

		PHASE,
		LATENT_HEAT,

		SIZE
	};

	static size_t offset;

	void computeInitialTemperature(const Step &step, std::vector<std::vector<double> > &data);

	double computeHTC(
			const ConvectionConfiguration &convection, const Element *e, const Point &p, Step step,
			double temp) const;

	void convectionMatParameters(
			const ConvectionConfiguration &convection, const Element *e, const Point &p, Step step,
			double temp, double T_EXT,
			double &rho, double &dynamic_viscosity, double &dynamic_viscosity_T, double &heat_capacity, double &thermal_conductivity) const;

	const HeatTransferConfiguration &_configuration;
	const ResultsSelectionConfiguration &_propertiesConfiguration;
};

}



#endif /* SRC_ASSEMBLER_PHYSICS_HEATTRANSFER_H_ */
