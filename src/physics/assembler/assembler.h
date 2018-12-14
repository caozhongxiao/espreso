
#ifndef SRC_PHYSICS_ASSEMBLER_ASSEMBLER_H_
#define SRC_PHYSICS_ASSEMBLER_ASSEMBLER_H_

namespace espreso {

struct HeatTransferLoadStepConfiguration;
struct FETISolverConfiguration;

struct Assembler {
	virtual void init() =0;
	virtual void nextTime() =0;
	virtual void assemble(Matrices matrices) =0;
	virtual void postProcess() =0;

	virtual ~Assembler() {};
};

template <typename TController, typename TComposer>
struct AssemblerInstance: public Assembler, public TController, public TComposer {

	AssemblerInstance(HeatTransferLoadStepConfiguration &loadStep, int DOFs)
	: TController(loadStep),
	  TComposer(*this, DOFs) {}

	AssemblerInstance(HeatTransferLoadStepConfiguration &loadStep, FETISolverConfiguration &solver, int DOFs)
	: TController(loadStep),
	  TComposer(*this, solver, DOFs) {}

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
		TComposer::assemble(matrices);
		TComposer::setDirichlet();
		TComposer::synchronize();
	}

	void postProcess()
	{
		TController::parametersChanged();
		TController::processSolution();
	}
};

}



#endif /* SRC_PHYSICS_ASSEMBLER_ASSEMBLER_H_ */
