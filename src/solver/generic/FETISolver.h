/*
 * LinearSolver.h
 *
 *  Created on: Sep 24, 2015
 *      Author: lriha
 */

#ifndef SOLVER_GENERIC_LINEARSOLVER_H_
#define SOLVER_GENERIC_LINEARSOLVER_H_

#include "linearsolver/linearsolver.h"
#include "solver/specific/itersolvers.h"
//#include "specific/superclusters.h"



namespace espreso {

struct DataHolder;

class FETISolver: public LinearSolver {
public:

	FETISolver(DataHolder *instance, FETISolverConfiguration &configuration);

	void init();

	void update(Matrices matrices);
	void solve();

	Matrices usedMatrices();
	double& precision() { return configuration.precision; }

	virtual ~FETISolver();

//	void setup();

	void init(const std::vector<int> &neighbours);

	void Preprocessing( std::vector < std::vector < esint > > & lambda_map_sub );

	void Solve( std::vector < std::vector <double > >	& f_vec, std::vector < std::vector < double > > & prim_solution );
	void Solve( std::vector < std::vector <double > >	& f_vec, std::vector < std::vector < double > > & prim_solution, std::vector < std::vector < double > > & dual_solution );

	void Postprocessing ();

	void finalize();

	void CheckSolution( std::vector < std::vector < double > > & prim_solution );

	void createCMat();

	DataHolder *instance;
	FETISolverConfiguration &configuration;

	TimeEval timeEvalMain;

private:

	SuperCluster *cluster;
	IterSolver   *solver;

	void setup_HTFETI();

	void setup_LocalSchurComplement();
	void setup_Preconditioner();
	void setup_FactorizationOfStiffnessMatrices();
	void setup_SetDirichletBoundaryConditions();

	void setup_CreateG_GGt_CompressG();
	void setup_InitClusterAndSolver();
};

}

#endif /* SOLVER_GENERIC_LINEARSOLVER_H_ */
