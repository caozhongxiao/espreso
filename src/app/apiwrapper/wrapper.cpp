
#include "wrapper.h"
#include "../../../include/feti4i.h"

#include "../../assembler/physics/precomputed.h"
#include "../../assembler/physicssolver/assembler.h"
#include "../../assembler/physicssolver/timestep/linear.h"
#include "../../assembler/physicssolver/loadstep/steadystate.h"
#include "../../assembler/step.h"

#include "../../config/ecf/ecf.h"
#include "../../input/api/api.h"

#include "../../mesh/mesh.h"
#include "../../output/result/resultstore.h"
#include "../../solver/generic/FETISolver.h"


espreso::ECFConfiguration* espreso::DataHolder::configuration = NULL;
std::list<FETI4IStructMatrix*> espreso::DataHolder::matrices;
std::list<FETI4IStructInstance*> espreso::DataHolder::instances;
espreso::TimeEval espreso::DataHolder::timeStatistics("API total time");

using namespace espreso;

FETI4IStructInstance::FETI4IStructInstance(FETI4IStructMatrix &matrix, eslocal *l2g, size_t size)
: instance(NULL), physics(NULL), linearSolver(NULL), assembler(NULL), timeStepSolver(NULL), loadStepSolver(NULL)
{
	store = new ResultStore();
	// TODO: MESH
	//mesh = new APIMesh(l2g, size);
}

FETI4IStructInstance::~FETI4IStructInstance()
{
	if (instance != NULL) { delete instance; }
	if (physics != NULL) { delete physics; }
	if (linearSolver != NULL) { delete linearSolver; }
	if (step != NULL) { delete step; }
	if (assembler != NULL) { delete assembler; }
	if (timeStepSolver != NULL) { delete timeStepSolver; }
	if (loadStepSolver != NULL) { delete loadStepSolver; }

	delete store;
	delete mesh;
}

static void checkConfiguration()
{
	if (espreso::DataHolder::configuration == NULL) {
		espreso::DataHolder::configuration = new ECFConfiguration();
		std::ifstream is("espreso.ecf");
		if (is.good()) {
			espreso::ECFReader::read(*espreso::DataHolder::configuration, "espreso.ecf", espreso::DataHolder::configuration->default_args, espreso::DataHolder::configuration->variables);
			espreso::ECFReader::set(espreso::DataHolder::configuration->environment, espreso::DataHolder::configuration->output);
		}
	}
}

void FETI4ISetDefaultIntegerOptions(FETI4IInt* options)
{
	checkConfiguration();
	ECFConfiguration &ecf = *espreso::DataHolder::configuration;

	options[FETI4I_SUBDOMAINS] = ecf.feti4ilibrary.domains;

	options[FETI4I_MAX_ITERATIONS] = ecf.feti4ilibrary.solver.max_iterations;
	options[FETI4I_FETI_METHOD] = static_cast<int>(ecf.feti4ilibrary.solver.method);
	options[FETI4I_PRECONDITIONER] = static_cast<int>(ecf.feti4ilibrary.solver.preconditioner);
	options[FETI4I_CGSOLVER] = static_cast<int>(ecf.feti4ilibrary.solver.iterative_solver);
	options[FETI4I_N_MICS] = ecf.feti4ilibrary.solver.n_mics;
	options[FETI4I_SC_SIZE] = ecf.feti4ilibrary.solver.sc_size;

	options[FETI4I_VERBOSE_LEVEL] = ecf.environment.verbose_level;
	options[FETI4I_MEASURE_LEVEL] = ecf.environment.measure_level;
	options[FETI4I_PRINT_MATRICES] = ecf.environment.print_matrices;
}

void FETI4ISetDefaultRealOptions(FETI4IReal* options)
{
	checkConfiguration();
	ECFConfiguration &ecf = *espreso::DataHolder::configuration;
	options[FETI4I_PRECISION] = ecf.feti4ilibrary.solver.precision;
}

static void FETI4ISetIntegerOptions(ECFConfiguration &configuration, FETI4IInt* options)
{
	if (!configuration.feti4ilibrary.getParameter(&configuration.feti4ilibrary.domains)->setValue(std::to_string(options[FETI4I_SUBDOMAINS]))) {
		ESINFO(GLOBAL_ERROR) << "Cannot set parameter 'SUBDOMAINS' to " << options[FETI4I_SUBDOMAINS];
	}
	if (!configuration.feti4ilibrary.solver.getParameter(&configuration.feti4ilibrary.solver.max_iterations)->setValue(std::to_string(options[FETI4I_MAX_ITERATIONS]))) {
		ESINFO(GLOBAL_ERROR) << "Cannot set parameter 'MAX_ITERATIONS' to " << options[FETI4I_MAX_ITERATIONS];
	}
	if (!configuration.feti4ilibrary.solver.getParameter(&configuration.feti4ilibrary.solver.method)->setValue(std::to_string(options[FETI4I_FETI_METHOD]))) {
		ESINFO(GLOBAL_ERROR) << "Cannot set parameter 'FETI_METHOD' to " << options[FETI4I_FETI_METHOD];
	}
	if (!configuration.feti4ilibrary.solver.getParameter(&configuration.feti4ilibrary.solver.preconditioner)->setValue(std::to_string(options[FETI4I_PRECONDITIONER]))) {
		ESINFO(GLOBAL_ERROR) << "Cannot set parameter 'PRECONDITIONER' to " << options[FETI4I_PRECONDITIONER];
	}
	if (!configuration.feti4ilibrary.solver.getParameter(&configuration.feti4ilibrary.solver.iterative_solver)->setValue(std::to_string(options[FETI4I_CGSOLVER]))) {
		ESINFO(GLOBAL_ERROR) << "Cannot set parameter 'CGSOLVER' to " << options[FETI4I_CGSOLVER];
	}
	if (!configuration.feti4ilibrary.solver.getParameter(&configuration.feti4ilibrary.solver.n_mics)->setValue(std::to_string(options[FETI4I_N_MICS]))) {
		ESINFO(GLOBAL_ERROR) << "Cannot set parameter 'n_mics' to " << options[FETI4I_N_MICS];
	}
	if (!configuration.feti4ilibrary.solver.getParameter(&configuration.feti4ilibrary.solver.sc_size)->setValue(std::to_string(options[FETI4I_SC_SIZE]))) {
		ESINFO(GLOBAL_ERROR) << "Cannot set parameter 'sc_size' to " << options[FETI4I_SC_SIZE];
	}
	if (!configuration.environment.getParameter(&configuration.environment.verbose_level)->setValue(std::to_string(options[FETI4I_VERBOSE_LEVEL]))) {
		ESINFO(GLOBAL_ERROR) << "Cannot set parameter 'VERBOSE_LEVEL' to " << options[FETI4I_VERBOSE_LEVEL];
	}
	if (!configuration.environment.getParameter(&configuration.environment.measure_level)->setValue(std::to_string(options[FETI4I_MEASURE_LEVEL]))) {
		ESINFO(GLOBAL_ERROR) << "Cannot set parameter 'MEASURE_LEVEL' to " << options[FETI4I_MEASURE_LEVEL];
	}
	if (!configuration.environment.getParameter(&configuration.environment.print_matrices)->setValue(std::to_string(options[FETI4I_PRINT_MATRICES]))) {
		ESINFO(GLOBAL_ERROR) << "Cannot set parameter 'PRINT_MATRICES' to " << options[FETI4I_PRINT_MATRICES];
	}

	configuration.feti4ilibrary.solver.regularization = FETI_REGULARIZATION::ALGEBRAIC;
}

static void FETI4ISetRealOptions(espreso::ECFConfiguration &configuration, FETI4IReal* options)
{
	configuration.feti4ilibrary.solver.precision = options[FETI4I_PRECISION];
}

void FETI4ICreateStiffnessMatrix(
		FETI4IMatrix 	*matrix,
		FETI4IInt		type,
		FETI4IInt		indexBase)
{
	checkConfiguration();
	DataHolder::timeStatistics.totalTime.startWithBarrier();
	TimeEvent event("Add element");
	DataHolder::timeStatistics.addEvent(event);

	DataHolder::matrices.push_back(new FETI4IStructMatrix(type, indexBase));
	*matrix = DataHolder::matrices.back();
}

void FETI4IAddElement(
		FETI4IMatrix	matrix,
		FETI4IInt		type,
		FETI4IInt		nodesSize,
		FETI4IInt*		nodes,
		FETI4IInt		dofsSize,
		FETI4IInt*		dofs,
		FETI4IReal*		values)
{
	espreso::DataHolder::timeStatistics.timeEvents.back().startWithoutBarrier();

	if (std::all_of(values, values + dofsSize, [] (double &value) { return value == 0; })) {
		// Skip elements with zero values
		return;
	}

	eslocal offset = matrix->offset;
	matrix->eType.push_back(type);
	matrix->eNodes.push_back(std::vector<eslocal>(nodes, nodes + nodesSize));
	matrix->eDOFs.push_back(std::vector<eslocal>(dofs, dofs + dofsSize));
	std::for_each(matrix->eNodes.back().begin(), matrix->eNodes.back().end(), [ &offset ] (eslocal &index) { index -= offset; });
	std::for_each(matrix->eDOFs.back().begin(), matrix->eDOFs.back().end(), [ &offset ] (eslocal &index) { index -= offset; });
	matrix->eMatrices.push_back(std::vector<double>(values, values + dofsSize * dofsSize));

	espreso::DataHolder::timeStatistics.timeEvents.back().endWithoutBarrier();
}

void FETI4ICreateInstance(
		FETI4IInstance 	*instance,
		FETI4IMatrix 	matrix,
		FETI4IInt 		size,
		FETI4IReal* 	rhs,
		FETI4IInt* 		l2g,
		FETI4IMPIInt 	neighbours_size,
		FETI4IMPIInt*	neighbours,
		FETI4IInt 		dirichlet_size,
		FETI4IInt* 		dirichlet_indices,
		FETI4IReal* 	dirichlet_values,
		FETI4IInt* 		integer_options,
		FETI4IReal*		real_options)
{
	checkConfiguration();
	DataHolder::instances.push_back(new FETI4IStructInstance(*matrix, l2g, size));

	FETI4ISetIntegerOptions(DataHolder::instances.back()->configuration, integer_options);
	FETI4ISetRealOptions(DataHolder::instances.back()->configuration, real_options);

	TimeEvent event("Create FETI4I instance"); event.startWithBarrier();

	ESINFO(OVERVIEW) << "ESPRESO create solver instance";

	std::vector<int> neighClusters = std::vector<int>(neighbours, neighbours + neighbours_size);

	input::API::load(
			*DataHolder::instances.back()->mesh, matrix->offset, DataHolder::instances.back()->configuration.feti4ilibrary.domains,
			matrix->eType, matrix->eNodes, matrix->eDOFs, matrix->eMatrices,
			dirichlet_size, dirichlet_indices, dirichlet_values,
			neighClusters,
			size, l2g);

	*instance = DataHolder::instances.back();

	DataHolder::instances.back()->instance = new Instance(*DataHolder::instances.back()->mesh);
	DataHolder::instances.back()->physics = new Precomputed(DataHolder::instances.back()->mesh, DataHolder::instances.back()->instance, (espreso::MatrixType)matrix->type, rhs, size);
	DataHolder::instances.back()->linearSolver = new FETISolver(DataHolder::instances.back()->instance, DataHolder::instances.back()->configuration.feti4ilibrary.solver);
	DataHolder::instances.back()->step = new Step();
	DataHolder::instances.back()->assembler = new Assembler(
			*DataHolder::instances.back()->instance,
			*DataHolder::instances.back()->physics,
			*DataHolder::instances.back()->mesh,
			*DataHolder::instances.back()->step,
			*DataHolder::instances.back()->store,
			*DataHolder::instances.back()->linearSolver);
	DataHolder::instances.back()->timeStepSolver = new LinearTimeStep(*DataHolder::instances.back()->assembler);
	DataHolder::instances.back()->loadStepSolver = new SteadyStateSolver(*DataHolder::instances.back()->timeStepSolver, 1);
	DataHolder::instances.back()->physics->prepare();

	event.endWithBarrier(); DataHolder::timeStatistics.addEvent(event);
}

void FETI4ISolve(
		FETI4IInstance 	instance,
		FETI4IInt 		solution_size,
		FETI4IReal*		solution)
{
	checkConfiguration();
	TimeEvent event("Solve FETI4I instance"); event.startWithBarrier();

	Logging::step = instance->step;
	instance->loadStepSolver->run();

	// TODO: MESH
	// memcpy(solution, instance->instance->solutions[espreso::Precomputed::SolutionIndex::MERGED]->data[0].data(), solution_size * sizeof(double));

	event.endWithBarrier(); DataHolder::timeStatistics.addEvent(event);
	DataHolder::timeStatistics.totalTime.endWithBarrier();
	DataHolder::timeStatistics.printStatsMPI();
}

template <typename TFETI4I>
static void destroy(std::list<TFETI4I*> &list, void *value)
{
	for (typename std::list<TFETI4I*>::iterator it = list.begin(); it != list.end(); ++it) {
		if (*it == value) {
			delete *it;
			list.erase(it);
			return;
		}
	}
}

void FETI4IDestroy(void *data)
{
	destroy(DataHolder::matrices, data);
	destroy(DataHolder::instances, data);
	if (DataHolder::instances.size() == 0 && DataHolder::matrices.size() == 0) {
		if (&DataHolder::configuration->environment == environment) {
			environment = NULL;
		}
		delete DataHolder::configuration;
		DataHolder::configuration = NULL;
	}
}

