
#ifndef ESPRESO_WRAPPER_H_
#define ESPRESO_WRAPPER_H_

#include <iostream>
#include <list>
#include <vector>

#include "basis/logging/timeeval.h"
#include "config/ecf/root.h"

namespace espreso {
struct Environment;
class DataHolder;
class Physics;
class TimeStepSolver;
class LoadStepSolver;
class Assembler;
class Mesh;
class ECFRoot;
class FETISolver;
class ResultStore;
}

struct FETI4IStructMatrix {
	FETI4IStructMatrix(esint type, esint offset): type(type), offset(offset) {};

	std::vector<esint> eType;
	std::vector<std::vector<esint> > eNodes;
	std::vector<std::vector<esint> > eDOFs;
	std::vector<std::vector<double> > eMatrices;

	esint type;
	esint offset;
};

struct FETI4IStructInstance {
	FETI4IStructInstance(FETI4IStructMatrix &matrix, esint *l2g, size_t size);
	~FETI4IStructInstance();

	espreso::DataHolder *instance;
	espreso::Physics * physics;
	espreso::FETISolver *linearSolver;
	espreso::ResultStore *store;
	espreso::Step *step;
	espreso::Assembler *assembler;
	espreso::TimeStepSolver *timeStepSolver;
	espreso::LoadStepSolver *loadStepSolver;

	espreso::Mesh *mesh;
	espreso::ECFRoot configuration;
};

namespace espreso {

struct APIDataHolder {
	static ECFRoot *configuration;
	static std::list<FETI4IStructMatrix*> matrices;
	static std::list<FETI4IStructInstance*> instances;
	static TimeEval timeStatistics;
};

}


#endif /* ESPRESO_WRAPPER_H_ */
