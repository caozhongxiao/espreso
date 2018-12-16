
#ifndef SRC_PHYSICS_ASSEMBLER_ASSEMBLER_H_
#define SRC_PHYSICS_ASSEMBLER_ASSEMBLER_H_

namespace espreso {

struct HeatTransferLoadStepConfiguration;
struct FETISolverConfiguration;

struct SolverParameters {
	double timeIntegrationConstantK;
	double timeIntegrationConstantM;
	bool tangentMatrixCorrection;
};

struct Assembler {
	virtual void init() =0;
	virtual void nextTime() =0;
	virtual void assemble(Matrices matrices, const SolverParameters &parameters) =0;
	virtual void postProcess() =0;

	virtual ~Assembler() {};
};

template <typename TController, typename TComposer, typename TProvider>
struct AssemblerInstance: public Assembler, public TController, public TComposer, public TProvider {

	AssemblerInstance(HeatTransferLoadStepConfiguration &loadStep)
	: TController(loadStep),
	  TComposer(*this, *this, 1),
	  TProvider(loadStep) {}

	AssemblerInstance(HeatTransferLoadStepConfiguration &loadStep, FETISolverConfiguration &solver)
	: TController(loadStep),
	  TComposer(*this, *this, solver, 1),
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

	void assemble(Matrices matrices, const SolverParameters &parameters)
	{
		TComposer::assemble(matrices, parameters);
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
