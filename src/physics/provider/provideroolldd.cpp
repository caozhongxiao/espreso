

#include "provideroolldd.h"

#include "mpi.h"

#include "../assembler/composer/composer.h"

#include "../../config/ecf/root.h"
#include "../../solver/generic/SparseMatrix.h"
#include "../../basis/logging/timeeval.h"
#include "../../basis/logging/logging.h"
#include "../../basis/utilities/utils.h"
#include "../../mesh/mesh.h"
#include "../../mesh/store/elementstore.h"

#include "../../output/result/resultstore.h"
#include "../../output/result/visualization/separated/vtklegacy.h"

#include "../../linearsolver/linearsolver.h"
#include "../dataholder.h"

using namespace espreso;

std::string ProviderOOLLDD::mNames(espreso::Matrices matrices, const std::string &prefix)
{
	return
	std::string(matrices & espreso::Matrices::K           ? prefix + "K "           : "") +
	std::string(matrices & espreso::Matrices::N           ? prefix + "N1 "          : "") +
	std::string(matrices & espreso::Matrices::N           ? prefix + "N2 "          : "") +
	std::string(matrices & espreso::Matrices::N           ? prefix + "RegMat "      : "") +
	std::string(matrices & espreso::Matrices::M           ? prefix + "M "           : "") +
	std::string(matrices & espreso::Matrices::R           ? prefix + "R "           : "") +
	std::string(matrices & espreso::Matrices::f           ? prefix + "f "           : "") +
	std::string(matrices & espreso::Matrices::B0          ? prefix + "B0 "          : "") +
	std::string(matrices & espreso::Matrices::B1          ? prefix + "B1 "          : "") +
	std::string(matrices & espreso::Matrices::B1c         ? prefix + "B1c "         : "") +
	std::string(matrices & espreso::Matrices::B1duplicity ? prefix + "B1duplicity " : "") +
	std::string(matrices & espreso::Matrices::primal      ? prefix + "Primal "      : "") +
	std::string(matrices & espreso::Matrices::dual        ? prefix + "Dual "        : "");
}

ProviderOOLLDD::ProviderOOLLDD(DataHolder &instance, Composer &composer, Mesh &mesh, LinearSolver &linearSolver)
: instance(instance), composer(composer), mesh(mesh), linearSolver(linearSolver), _timeStatistics(new TimeEval("Physics solver timing"))
{
	_timeStatistics->totalTime.startWithBarrier();

	composer.initDOFs();
	composer.initDirichlet();
	composer.buildPatterns();
}

void ProviderOOLLDD::preprocessData()
{
	timeWrapper("pre-process data", [&] () {
		composer.initData();
	});
}

void ProviderOOLLDD::updateStructuralMatrices(Matrices matrices)
{
	Matrices updated = matrices & (Matrices::K | Matrices::M | Matrices::f | Matrices::R);

	if (updated) {
		timeWrapper("update " + mNames(updated), [&] () {
//			physics.updateMatrix(updated);
//			composer.assemble(updated);
		});
	}
}

void ProviderOOLLDD::updateGluingMatrices(Matrices matrices)
{
	timeWrapper("update " + mNames(Matrices::B1c), [&] () {
		composer.setDirichlet();
	});

	timeWrapper("update " + mNames(Matrices::B1duplicity), [&] () {
		composer.synchronize();
	});
}

void ProviderOOLLDD::solve(Matrices updatedMatrices)
{
	Matrices solverMatrices = Matrices::K | Matrices::M | Matrices::f | Matrices::B1;
	storeWrapper(mNames(solverMatrices), solverMatrices);

	timeWrapper("update linear solver: " + mNames(updatedMatrices), [&] () {
		linearSolver.update(updatedMatrices);
	});

	timeWrapper("run linear solver", [&] () {
		linearSolver.solve();
	});

	timeWrapper("fill solution", [&] () {
		composer.fillSolution();
	});
}

void ProviderOOLLDD::nextTime()
{
	timeWrapper("update time", [&] () {
		composer.nextTime();
	});
}

void ProviderOOLLDD::parametersChanged()
{
	timeWrapper("update parameters", [&] () {
		composer.parametersChanged();
	});
}

void ProviderOOLLDD::processSolution()
{
	timeWrapper("post-processing", [&] () {
		composer.processSolution();
	});
	storeWrapper(mNames(Matrices::primal), Matrices::primal);
}

void ProviderOOLLDD::storeSolution()
{
	if (mesh.store->storeStep()) {
		timeWrapper("store solution", [&] () {
			mesh.store->updateSolution();
		});
	}
}

ProviderOOLLDD::~ProviderOOLLDD()
{
	delete _timeStatistics;
	for (auto it = _timeEvents.begin(); it != _timeEvents.end(); ++it) {
		delete it->second;
	}
}

void ProviderOOLLDD::finalize()
{
	timeWrapper("finalize", [&] () {
		linearSolver.finalize();
	});
	_timeStatistics->totalTime.endWithBarrier();
	_timeStatistics->printStatsMPI();
}

void ProviderOOLLDD::timeWrapper(const std::string &action, std::function<void(void)> operations)
{
	std::string fulldesc(/*physics.name() + ": " +*/ action);

	ESINFO(PROGRESS2) << fulldesc;

	TimeEvent *event;
	if (_timeEvents.find(fulldesc) != _timeEvents.end()) {
		event = _timeEvents[fulldesc];
	} else {
		_timeEvents[fulldesc] = event = new TimeEvent(fulldesc);
		_timeStatistics->addPointerToEvent(event);
	}

	event->start();
	operations();
	event->endWithBarrier();
}

template<typename TType>
static void storeData(TType &data, size_t domain, const std::string &name)
{
	std::ofstream os(Logging::prepareFile(domain, name));
	os.precision(10);
	os << data;
	os.close();
}


bool ProviderOOLLDD::checkForStore(const std::string &name)
{
	if (environment->print_matrices) {
		std::string fulldesc(/*physics.name() + ": store " +*/ name);
		ESINFO(ALWAYS_ON_ROOT) << Info::TextColor::BLUE << fulldesc;
	}
	return environment->print_matrices;
}

void ProviderOOLLDD::storeMatrices(Matrices matrices, size_t domain)
{
	auto storeMatrix = [&] (std::vector<SparseMatrix> &data, Matrices matrix, const std::string &name) {
		if (matrices & matrix) {
			storeData(data[domain], domain, name);
		}
	};

	auto storeVector = [&] (std::vector<std::vector<double> > &data, Matrices matrix, const std::string &name) {
		if (matrices & matrix) {
			storeData(data[domain], domain, name);
		}
	};

	storeMatrix(instance.K, Matrices::K, "K");
	storeMatrix(instance.N1, Matrices::N, "N1");
	storeMatrix(instance.N2, Matrices::N, "N2");
	storeMatrix(instance.RegMat, Matrices::N, "RegMat");
	storeMatrix(instance.M, Matrices::M, "M");
	storeVector(instance.R, Matrices::R, "R");
	storeVector(instance.f, Matrices::f, "f");
	storeMatrix(instance.B0, Matrices::B0, "B0");
	storeMatrix(instance.B1, Matrices::B1, "B1");
	storeVector(instance.B1c, (Matrices::B1 | Matrices::B1c), "B1c");
	storeVector(instance.B1duplicity, (Matrices::B1 | Matrices::B1duplicity), "B1duplicity");
	storeVector(instance.primalSolution, Matrices::primal, "solution");
	storeVector(instance.dualSolution, Matrices::dual, "dualSolution");
}

void ProviderOOLLDD::storeWrapper(const std::string &name, Matrices matrices)
{
	if (checkForStore(name)) {
		for (size_t d = 0; d < mesh.elements->ndomains; d++) {
			storeMatrices(matrices, d);
		}
	}
}

void ProviderOOLLDD::storeWrapper(const std::string &name, Matrices matrices, size_t domain)
{
	if (checkForStore(name)) {
		storeMatrices(matrices, domain);
	}
}

void ProviderOOLLDD::storeWrapper(const std::string &name, std::vector<SparseMatrix> &matrices)
{
	if (checkForStore(name)) {
		for (size_t d = 0; d < matrices.size(); d++) {
			storeData(matrices[d], d, name);
		}
	}
}

void ProviderOOLLDD::storeWrapper(const std::string &name, std::vector<std::vector<double> > &data)
{
	if (checkForStore(name)) {
		for (size_t d = 0; d < data.size(); d++) {
			storeData(data[d], d, name);
		}
	}
}








