
#ifndef SRC_PHYSICS_ASSEMBLER_ASSEMBLER_H_
#define SRC_PHYSICS_ASSEMBLER_ASSEMBLER_H_

#include <string>
#include <vector>

namespace espreso {

struct HeatTransferLoadStepConfiguration;
struct FETISolverConfiguration;
enum Matrices: int;
enum class SumRestriction;
class SparseMatrix;

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
	virtual void postProcess() =0;

	/// z = a * x + b + y
	void sum(std::vector<std::vector<double> > &z, double a, const std::vector<std::vector<double> > &x, double b, const std::vector<std::vector<double> > &y, const std::string &description);
	/// z = a * x + b + y (prefix variant)
	void sum(std::vector<std::vector<double> > &z, double a, const std::vector<std::vector<double> > &x, double b, const std::vector<std::vector<double> > &y, const std::vector<size_t> &prefix, const std::string &description);
	/// A += beta * B
	void sum(std::vector<SparseMatrix> &A, double beta, std::vector<SparseMatrix> &B, const std::string &description);

	/// y = A * x
	void multiply(std::vector<std::vector<double> > &y, std::vector<SparseMatrix> &A, std::vector<std::vector<double> > &x, const std::string &description);
	/// a = x * y
	double multiply(std::vector<std::vector<double> > &x, std::vector<std::vector<double> > &y, const std::string &description);

	double sumSquares(const std::vector<std::vector<double> > &data, SumRestriction restriction, const std::string &description);
	void addToDirichletInB1(double a, const std::vector<std::vector<double> > &x);
	double maxAbsValue(const std::vector<std::vector<double> > &v, const std::string &description);
	double lineSearch(const std::vector<std::vector<double> > &U, std::vector<std::vector<double> > &deltaU, std::vector<std::vector<double> > &F_ext);
	void keepK();

	virtual ~Assembler() {};

	SolverParameters parameters;
};

template <typename TController, typename TComposer, typename TProvider>
struct AssemblerInstance: public Assembler, public TController, public TComposer, public TProvider {

	template <typename TPhysics>
	AssemblerInstance(TPhysics &loadStep, int DOFs)
	: TController(loadStep),
	  TComposer(*this, *this, DOFs),
	  TProvider(loadStep) {}

	template <typename TPhysics>
	AssemblerInstance(TPhysics &loadStep, FETISolverConfiguration &solver, int DOFs)
	: TController(loadStep),
	  TComposer(*this, *this, solver, DOFs),
	  TProvider(loadStep) {}

	void init()
	{
		TComposer::initDOFs();
		TComposer::initDirichlet();
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
			keepK();
		}
	}

	void setDirichlet()
	{
		TComposer::setDirichlet();
		TComposer::synchronize();
	}

	void postProcess()
	{
		TComposer::fillSolution();

		TController::parametersChanged();
		TController::processSolution();
	}
};

}



#endif /* SRC_PHYSICS_ASSEMBLER_ASSEMBLER_H_ */
