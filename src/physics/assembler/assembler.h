
#ifndef SRC_PHYSICS_ASSEMBLER_ASSEMBLER_H_
#define SRC_PHYSICS_ASSEMBLER_ASSEMBLER_H_

#include <string>
#include <vector>

namespace espreso {

struct HeatTransferLoadStepConfiguration;
struct FETISolverConfiguration;
struct MultigridConfiguration;
enum Matrices: int;
class Composer;
class Controller;
class Provider;
struct NodeData;
enum class SumRestriction;
class SparseMatrix;
class LinearSolver;

struct SolverParameters {
	double internalForceReduction;
	double timeIntegrationConstantK;
	double timeIntegrationConstantM;
	bool tangentMatrixCorrection;

	SolverParameters()
	: internalForceReduction(1),
	  timeIntegrationConstantK(1),
	  timeIntegrationConstantM(0),
	  tangentMatrixCorrection(false)
	{}
};

struct Assembler {
	virtual void init() =0;
	virtual void nextTime() =0;
	virtual void parametersChanged() =0;
	virtual void assemble(Matrices matrices) =0;
	virtual void setDirichlet(Matrices matrices, const std::vector<double> &subtraction = {}) =0;
	virtual void solve(Matrices matrices) =0;
	virtual void postProcess() =0;

	virtual double& solutionPrecision() =0;

	double lineSearch(const std::vector<std::vector<double> > &U, std::vector<std::vector<double> > &deltaU, std::vector<std::vector<double> > &F_ext);

	virtual DataHolder* data() =0;
	virtual Composer* composer() =0;
	virtual Controller* controller() =0;
	virtual Provider* provider() =0;

	virtual ~Assembler() {};

	SolverParameters parameters;
	LinearSolver *_solver;

protected:
	void callsolve(Matrices matrices);
};

template <typename TController, typename TComposer, typename TProvider, typename TSolver>
struct AssemblerInstance: public Assembler, public TController, public TComposer, public TProvider {

	template <typename TPhysics>
	AssemblerInstance(TPhysics &loadStep, MultigridConfiguration &solver, int DOFs)
	: TController(loadStep),
	  TComposer(*this, *this, DOFs),
	  TProvider(TComposer::data, loadStep)
	{
		_solver = new TSolver(TComposer::data, solver);
	}

	template <typename TPhysics>
	AssemblerInstance(TPhysics &loadStep, FETISolverConfiguration &solver, int DOFs)
	: TController(loadStep),
	  TComposer(*this, *this, solver, DOFs),
	  TProvider(TComposer::data, loadStep)
	{
		_solver = new TSolver(TComposer::data, solver);
	}

	~AssemblerInstance()
	{
		if (_solver != NULL) {
			delete _solver;
		}
	}

	DataHolder* data() { return TComposer::data; }
	Composer* composer() { return this; }
	Controller* controller() { return this; }
	Provider* provider() { return this; }

	void init()
	{
		TComposer::initDOFs();
		TComposer::buildDirichlet();
		TComposer::buildPatterns();
		if (TProvider::needMatrixVectorProduct()) {
			TComposer::buildMVData();
		}

		TController::initData();
	}

	void nextTime()
	{
		TController::nextTime();
	}

	void parametersChanged()
	{
		TController::parametersChanged();
	}

	void assemble(Matrices matrices)
	{
		TComposer::assemble(matrices, parameters);
		if ((matrices & Matrices::K) && TProvider::needOriginalStiffnessMatrices()) {
			TComposer::keepK();
		}
		if ((matrices & Matrices::f) && TProvider::needOriginalRHS()) {
			TComposer::keepRHS();
		}
	}

	void setDirichlet(Matrices matrices, const std::vector<double> &subtraction = {})
	{
		TComposer::setDirichlet(matrices, subtraction);
	}

	void solve(Matrices matrices)
	{
		callsolve(matrices);
		TComposer::fillSolution();
		if (TProvider::needReactionForces()) {
			TComposer::computeReactionForces();
		}
	}

	void postProcess()
	{
		TComposer::processSolution();
	}

	double& solutionPrecision()
	{
		return TProvider::solutionPrecision();
	}
};

}



#endif /* SRC_PHYSICS_ASSEMBLER_ASSEMBLER_H_ */
