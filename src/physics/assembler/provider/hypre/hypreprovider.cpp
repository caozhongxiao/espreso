
#include "hypreprovider.h"

#include "basis/matrices/matrixtype.h"
#include "esinfo/eslog.hpp"
#include "config/ecf/physics/physicssolver/loadstep.h"

using namespace espreso;

HYPREProvider::HYPREProvider(DataHolder *data, LoadStepConfiguration &configuration)
: Provider(data, configuration)
{

}

MatrixType HYPREProvider::getMatrixType() const
{
	return MatrixType::REAL_UNSYMMETRIC; // Hypre always get full matrix
}

bool HYPREProvider::needMatrixVectorProduct()
{
	return Provider::needMatrixVectorProduct() || _configuration.mode == LoadStepConfiguration::MODE::NONLINEAR;
}

bool HYPREProvider::needOriginalStiffnessMatrices()
{
	return true;
//	return Provider::needOriginalStiffnessMatrices() || _configuration.mode == LoadStepConfiguration::MODE::NONLINEAR;
}

double& HYPREProvider::solutionPrecision()
{
	switch (_configuration.hypre.solver_type) {
	case HypreConfiguration::SOLVER_TYPE::BiCGSTAB:
		return _configuration.hypre.bicgstab.relative_conv_tol;
	case HypreConfiguration::SOLVER_TYPE::BoomerAMG:
		return _configuration.hypre.boomeramg.convergence_tolerance;
	case HypreConfiguration::SOLVER_TYPE::CGNR:
		return _configuration.hypre.cgnr.relative_conv_tol;
	case HypreConfiguration::SOLVER_TYPE::FlexGMRES:
		return _configuration.hypre.flexgmres.relative_conv_tol;
	case HypreConfiguration::SOLVER_TYPE::GMRES:
		return _configuration.hypre.gmres.relative_conv_tol;
	case HypreConfiguration::SOLVER_TYPE::LGMRES:
		return _configuration.hypre.lgmres.relative_conv_tol;
	case HypreConfiguration::SOLVER_TYPE::PCG:
		return _configuration.hypre.pcg.relative_conv_tol;
	default:
		eslog::globalerror("Required precision of unknown solver.\n");
		exit(0);
	}
}





