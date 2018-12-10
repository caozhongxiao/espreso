
#include "loadstepiterator.h"
#include "dataholder.h"

#include "assembler/assembler.h"

#include "../basis/logging/logging.h"
#include "../config/ecf/physics/heattransfer.h"
#include "../config/ecf/physics/structuralmechanics.h"
#include "../globals/time.h"

#include "../mesh/mesh.h"

#include "../linearsolver/multigrid/multigrid.h"
#include "../solver/generic/FETISolver.h"

using namespace espreso;

DataHolder* LoadStepIterator::getDataHolder(Mesh &mesh)
{
	return new DataHolder(mesh);
}

Composer* LoadStepIterator::getComposer(HeatTransferLoadStepConfiguration &configuration)
{
//	switch () {
//	case DIMENSION::D2:
//		switch (_loadStep->solver) {
//		case LoadStepConfiguration::SOLVER::FETI:
//			_composer.push_back(new DomainsHeatTransfer2D(
//							*mesh, *_dataHolder, configuration, configuration, propertiesConfiguration));
//			break;
//		case LoadStepConfiguration::SOLVER::MULTIGRID:
//			_composer.push_back(new GlobalHeatTransfer2D(
//							*mesh, *_dataHolder, configuration, configuration, propertiesConfiguration));
//			break;
//		}
//
//		break;
//	case DIMENSION::D3:
//		switch (_loadStep->solver) {
//		case LoadStepConfiguration::SOLVER::FETI:
//			_composer.push_back(new DomainsHeatTransfer3D(
//							*mesh, *_instances.front(), *step, configuration, configuration.load_steps_settings.at(1), propertiesConfiguration));
//			break;
//		case LoadStepConfiguration::SOLVER::MULTIGRID:
//			_composer.push_back(new GlobalHeatTransfer3D(
//							*mesh, *_instances.front(), *step, configuration, configuration.load_steps_settings.at(1), propertiesConfiguration));
//			break;
//		}
//		break;
//	default:
//		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: invalid dimension.";
//	}
}

LinearSolver* LoadStepIterator::getLinearSolver()
{
	switch (_loadStep->solver) {
	case LoadStepConfiguration::SOLVER::FETI:
		return new FETISolver(_dataHolder, _loadStep->feti);
	case LoadStepConfiguration::SOLVER::MULTIGRID:
		return new MultigridSolver(_dataHolder, _loadStep->multigrid);
	default:
		ESINFO(GLOBAL_ERROR) << "Not implemented requested SOLVER.";
		return NULL;
	}
}

LoadStepIterator::LoadStepIterator(Mesh &mesh, HeatTransferConfiguration &configuration)
{
	_loadStep = &configuration.load_steps_settings.at(time::step + 1);

	_dataHolder = getDataHolder(mesh);
	_linearSolver = getLinearSolver();
//	_composer = getComposer();
}

//LoadStepIterator::LoadStepIterator(Mesh &mesh, StructuralMechanicsConfiguration &configuration)
//{
//
//}


