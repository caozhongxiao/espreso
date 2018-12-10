
#ifndef SRC_PHYSICS_LOADSTEPITERATOR_H_
#define SRC_PHYSICS_LOADSTEPITERATOR_H_

namespace espreso {

struct HeatTransferConfiguration;
struct StructuralMechanicsConfiguration;
struct LoadStepConfiguration;
struct HeatTransferLoadStepConfiguration;

class Mesh;

struct DataHolder;
class LoadStepSolver;
class TimeStepSolver;
class Composer;
class LinearSolver;


class LoadStepIterator {

public:
	LoadStepIterator(Mesh &mesh, HeatTransferConfiguration &configuration);
//	LoadStepIterator(Mesh &mesh, StructuralMechanicsConfiguration &configuration);

protected:
	DataHolder *_dataHolder;
	LoadStepSolver *_loadStepSolver;
	TimeStepSolver *_timeStepSolver;
	Composer *_composer;
	LinearSolver * _linearSolver;

	LoadStepConfiguration *_loadStep;


	DataHolder* getDataHolder(Mesh &mesh);
	Composer* getComposer(HeatTransferLoadStepConfiguration &configuration);
	LinearSolver* getLinearSolver();
};

}



#endif /* SRC_PHYSICS_LOADSTEPITERATOR_H_ */
