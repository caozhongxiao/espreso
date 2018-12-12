
#include "loadstepiterator.h"
#include "dataholder.h"

#include "assembler/assembler.h"

#include "../globals/run.h"
#include "../globals/time.h"
#include "../basis/logging/logging.h"
#include "../config/ecf/root.h"

#include "../mesh/mesh.h"

#include "../linearsolver/multigrid/multigrid.h"
#include "../solver/generic/FETISolver.h"

using namespace espreso;


//Composer* LoadStepIterator::getComposer(HeatTransferLoadStepConfiguration &configuration)
//{
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
//}

static LinearSolver* getLinearSolver(LoadStepConfiguration &loadStep)
{
	switch (loadStep.solver) {
	case LoadStepConfiguration::SOLVER::FETI:
		return new FETISolver(run::data, loadStep.feti);
	case LoadStepConfiguration::SOLVER::MULTIGRID:
		return new MultigridSolver(loadStep.multigrid);
	default:
		ESINFO(GLOBAL_ERROR) << "Not implemented requested SOLVER.";
		return NULL;
	}
}

//static Assembler* getAssembler(LoadStepConfiguration &loadStep)
//{
//	switch (loadStep.solver) {
//	case LoadStepConfiguration::SOLVER::FETI:
//		return new FETISolver(run::data, loadStep.feti);
//	case LoadStepConfiguration::SOLVER::MULTIGRID:
//		return new MultigridSolver(loadStep.multigrid);
//	default:
//		ESINFO(GLOBAL_ERROR) << "Not implemented requested SOLVER.";
//		return NULL;
//	}
//}

LoadStepIterator::LoadStepIterator()
: _loadStepSolver(NULL), _timeStepSolver(NULL), _assembler(NULL), _linearSolver(NULL)
{
	run::data = new DataHolder();
}

bool LoadStepIterator::next()
{
	switch (run::ecf->physics) {
	case PHYSICS::HEAT_TRANSFER_2D:
		next(run::ecf->heat_transfer_2d);
		break;
	case PHYSICS::HEAT_TRANSFER_3D:
		next(run::ecf->heat_transfer_3d);
		break;
	case PHYSICS::STRUCTURAL_MECHANICS_2D:
	case PHYSICS::STRUCTURAL_MECHANICS_3D:
		next();
	default:
		ESINFO(GLOBAL_ERROR) << "Unknown physics.";
	}

	return false;
}

bool LoadStepIterator::next(HeatTransferConfiguration &configuration)
{
	if (time::step++ < configuration.load_steps) {
		_linearSolver = getLinearSolver(configuration.load_steps_settings.at(time::step));
//		_assembler = getAssembler(configuration.load_steps_settings.at(time::step));
		_assembler = new GlobalAssembler<UniformNodesComposer, HeatTransfer2DControler>();
	}

	return time::step < configuration.load_steps;
}

bool LoadStepIterator::next(StructuralMechanicsConfiguration &configuration)
{
	if (time::step < configuration.load_steps) {



		++time::step;
	}

	return time::step < configuration.load_steps;
}


