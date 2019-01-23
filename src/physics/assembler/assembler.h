
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
	virtual void assemble(Matrices matrices) =0;
	virtual void setDirichlet() =0;
	virtual void solve(Matrices matrices) =0;
	virtual void postProcess() =0;

	virtual NodeData* solution() =0;
	virtual double& solutionPrecision() =0;

	/// z = a * x + b + y
	void sum(std::vector<std::vector<double> > &z, double a, const std::vector<std::vector<double> > &x, double b, const std::vector<std::vector<double> > &y, const std::string &description);
	/// z = a * x + b + y (prefix variant)
	void sum(std::vector<std::vector<double> > &z, double a, const std::vector<std::vector<double> > &x, double b, const std::vector<std::vector<double> > &y, const std::vector<size_t> &prefix, const std::string &description);

	/// y = A * x
	void multiply(std::vector<std::vector<double> > &y, std::vector<SparseMatrix> &A, std::vector<std::vector<double> > &x, const std::string &description);
	/// a = x * y
	double multiply(std::vector<std::vector<double> > &x, std::vector<std::vector<double> > &y, const std::string &description);

	double sumSquares(const std::vector<std::vector<double> > &data, SumRestriction restriction, const std::string &description);
	void addToDirichletInB1(double a, const std::vector<std::vector<double> > &x);
	double maxAbsValue(const std::vector<std::vector<double> > &v, const std::string &description);
	double lineSearch(const std::vector<std::vector<double> > &U, std::vector<std::vector<double> > &deltaU, std::vector<std::vector<double> > &F_ext);

	virtual DataHolder* data() =0;
	virtual Composer* composer() =0;

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

	void init()
	{
		TComposer::initDOFs();
		TComposer::buildDirichlet();
		TComposer::buildPatterns();

		TController::initData();
	}

	void nextTime()
	{
		TController::nextTime();
	}

	void assemble(Matrices matrices)
	{
		TComposer::assemble(matrices, parameters);
		if (TProvider::needOriginalStiffnessMatrices()) {
			TComposer::keepK();
		}
	}

	void setDirichlet()
	{
		TComposer::setDirichlet();
	}

	void solve(Matrices matrices)
	{
		callsolve(matrices);
		TComposer::fillSolution();
	}

	void postProcess()
	{
		TComposer::parametersChanged();
		TComposer::processSolution();
	}

	NodeData* solution()
	{
		return TController::solution();
	}

	double& solutionPrecision()
	{
		return TProvider::solutionPrecision();
	}
};

}



#endif /* SRC_PHYSICS_ASSEMBLER_ASSEMBLER_H_ */
