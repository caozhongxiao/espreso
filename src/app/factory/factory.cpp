
#include <signal.h>
#include <csignal>

#include "factory.h"

#include "heattransferfactory.h"
#include "structuralmechanicsfactory.h"

#include "../../assembler/physicssolver/timestep/timestepsolver.h"
#include "../../assembler/physicssolver/loadstep/loadstepsolver.h"
#include "../../assembler/physicssolver/assembler.h"
#include "../../assembler/physics/physics.h"
#include "../../assembler/step.h"
#include "../../assembler/instance.h"

#include "../../newmesh/newmesh.h"
#include "../../input/loader.h"
#include "../../config/ecf/ecf.h"
#include "../../mesh/structures/mesh.h"
#include "../../solver/generic/FETISolver.h"

#include "../../output/resultstorelist.h"

namespace espreso {

static void signalHandler(int signal)
{
	switch (signal) {
	case SIGSEGV:
		ESINFO(ERROR) << "Invalid memory reference";
		break;
	case SIGFPE:
		ESINFO(ERROR) << "Erroneous arithmetic operation";
		break;
	default:
		ESINFO(ERROR) << "ESPRESO trigger error " << signal << ".";
	}
}

void ESPRESO::run(int *argc, char ***argv)
{
	std::signal(SIGFPE, signalHandler);
	std::signal(SIGSEGV, signalHandler);

	ECFConfiguration configuration(argc, argv);
	ESINFO(OVERVIEW) << "Run ESPRESO on " << environment->MPIsize << " process(es).";

	Mesh mesh;
	ResultStoreList* solutionStore = ResultStoreList::createAsynchronizedStore(configuration.output);
	if (ResultStoreList::isComputeNode()) {
		Factory factory(configuration, mesh, *solutionStore);
		factory.solve();
	}
	ResultStoreList::destroyAsynchronizedStore();
}

Factory::Factory(const ECFConfiguration &configuration, Mesh &mesh, ResultStoreList &store)
: _mesh(&mesh), _store(&store), _loader(NULL)
{
	NewMesh nmesh(mesh);
	Loader::load(configuration, nmesh, configuration.environment.MPIrank, configuration.environment.MPIsize);

	// LOAD PHYSICS
	switch (configuration.physics) {
	case PHYSICS::HEAT_TRANSFER_2D:
		_loader = new HeatTransferFactory(configuration.heat_transfer_2d, configuration.output.results_selection, _mesh);
		break;
	case PHYSICS::HEAT_TRANSFER_3D:
		_loader = new HeatTransferFactory(configuration.heat_transfer_3d, configuration.output.results_selection, _mesh);
		break;
	case PHYSICS::STRUCTURAL_MECHANICS_2D:
		_loader = new StructuralMechanicsFactory(configuration.structural_mechanics_2d, configuration.output.results_selection, _mesh);
		break;
	case PHYSICS::STRUCTURAL_MECHANICS_3D:
		_loader = new StructuralMechanicsFactory(configuration.structural_mechanics_3d, configuration.output.results_selection, _mesh);
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Unknown PHYSICS in configuration file";
	}

	for (size_t step = 0; step < _loader->loadSteps(); step++) {
		_loadSteps.push_back(_loader->getLoadStepSolver(step, _mesh, _store));
	}

	_loader->preprocessMesh();
	_store->updateMesh();
}

void Factory::solve()
{
	Step step;
	Logging::step = &step;

	for (step.step = 0; step.step < _loadSteps.size(); step.step++) {
		_loadSteps[step.step]->run(step);
	}
}

template <class TType>
static void clear(std::vector<TType> &vector)
{
	for (size_t i = 0; i < vector.size(); i++) {
		delete vector[i];
	}
}

FactoryLoader::~FactoryLoader()
{
	clear(_instances);
	clear(_physics);
	clear(_linearSolvers);
	clear(_assemblers);
	clear(_timeStepSolvers);
	clear(_loadStepSolvers);
}

LinearSolver* FactoryLoader::getLinearSolver(const LoadStepConfiguration &settings, Instance *instance) const
{
	switch (settings.solver) {
	case LoadStepConfiguration::SOLVER::FETI:
		return new FETISolver(instance, settings.feti);
	default:
		ESINFO(GLOBAL_ERROR) << "Not implemented requested SOLVER.";
		return NULL;
	}
}

void FactoryLoader::preprocessMesh()
{
	// TODO: generalize it !!

	for (size_t i = 0; i < _physics.size(); i++) {

		switch (dynamic_cast<FETISolver*>(_linearSolvers.front())->configuration.method) {
		case FETI_METHOD::TOTAL_FETI:
			_physics[i]->prepare();
			break;
		case FETI_METHOD::HYBRID_FETI:
			switch (dynamic_cast<FETISolver*>(_linearSolvers.front())->configuration.B0_type) {
			case FETI_B0_TYPE::CORNERS:
				_physics[i]->prepareHybridTotalFETIWithCorners();
				break;
			case FETI_B0_TYPE::KERNELS:
				_physics[i]->prepareHybridTotalFETIWithKernels();
				break;
			default:
				ESINFO(GLOBAL_ERROR) << "Unknown type of B0";
			}
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "Unknown FETI method";
		}
	}
}

Factory::~Factory()
{
	delete _loader;
}

void FactoryLoader::printError(const std::string &error) const
{
	ESINFO(GLOBAL_ERROR) << error;
}

}


