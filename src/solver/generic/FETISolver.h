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

#include "physics/dataholder.h"


namespace espreso {

struct DataHolder;

class FETISolver: public LinearSolver {
public:

	FETISolver(DataHolder *instance, const FETISolverConfiguration &configuration);

	void init();

	void update(Matrices matrices);
	void solve();

	bool glueDomainsByLagrangeMultipliers() const { return true; }
	bool applyB1Scaling() const { return configuration.scaling; }
	bool applyB1LagrangeRedundancy() const { return configuration.redundant_lagrange; }

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
	FETISolverConfiguration configuration;

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
