
#include "../../include/feti4i.h"

#include "factory/factory.h"

#include "../config/ecf/ecf.h"
#include "../basis/logging/logging.h"
#include "../basis/matrices/denseMatrix.h"
#include "../mesh/elements/element.h"
#include "../mesh/structures/mesh.h"
#include "../mesh/structures/coordinates.h"
#include "../assembler/step.h"
#include "../assembler/instance.h"
#include "../assembler/constraints/equalityconstraints.h"
#include "../assembler/physics/physics.h"
#include "../assembler/physicssolver/assembler.h"
#include "../solver/generic/SparseMatrix.h"

namespace espreso {

class APITestESPRESODataProvider {

public:
	APITestESPRESODataProvider(int *argc, char ***argv): ecf(argc, argv), factory(ecf, 1)
	{
		if (factory._loadSteps.size() > 1) {
			ESINFO(GLOBAL_ERROR) << "APITEST: Cannot test instance with more loadsteps.";
		}
		factory._loader->_physics[0]->preprocessData(step);
	}

	int matrixType()
	{
		return static_cast<int>(factory._loader->_assemblers[0]->physics.getMatrixType(step, 0));
	}

	size_t elements()
	{
		return factory._mesh->elements().size();
	}

	void addElementMatrix(FETI4IMatrix &K, size_t e)
	{
		FETI4IInt eType = static_cast<int>(factory._mesh->elements()[e]->type());
		std::vector<FETI4IInt> nodes;
		for (size_t n = 0; n < factory._mesh->elements()[e]->nodes(); n++) {
			nodes.push_back(factory._mesh->elements()[e]->node(n));
		}
		std::vector<FETI4IInt> dofs;
		factory._loader->_physics[0]->fillDOFsIndices(factory._mesh->elements()[e], 0, dofs);
		DenseMatrix eK, eM, eR, eF;
		factory._loader->_physics[0]->updateMatrix(step, Matrices::K, factory._mesh->elements()[e], eK, eM, eR, eF, factory._loader->_instances[0]->solutions);
		FETI4IAddElement(K, eType, nodes.size(), nodes.data(), dofs.size(), dofs.data(), eK.values());
	}

	void computeRHS(std::vector<FETI4IReal> &rhs)
	{
		rhs.resize(factory._loader->_instances[0]->domainDOFCount[0]);
		std::vector<FETI4IInt> dofs;
		for (size_t e = 0; e < elements(); e++) {
			factory._loader->_physics[0]->fillDOFsIndices(factory._mesh->elements()[e], 0, dofs);
			DenseMatrix eK, eM, eR, eF;
			factory._loader->_physics[0]->updateMatrix(step, Matrices::f, factory._mesh->elements()[e], eK, eM, eR, eF, factory._loader->_instances[0]->solutions);
			for (size_t dof = 0; dof < dofs.size(); dof++) {
				rhs[dofs[dof]] = eF(0, dof);
			}
		}
	}

	void fillL2G(std::vector<FETI4IInt> &l2g)
	{
		for (size_t n = 0; n < factory._mesh->nodes().size(); n++) {
			if (factory._mesh->nodes()[n]->parentElements().size()) {
				l2g.push_back(factory._mesh->coordinates().globalIndex(n));
			}
		}
	}

	void fillDirichlet(std::vector<FETI4IInt> &dirichlet_indices, std::vector<FETI4IReal> &dirichlet_values)
	{
		factory._loader->_physics[0]->_equalityConstraints->insertDirichletToB1(step, false);
		for (auto it = factory._loader->_instances[0]->B1[0].J_col_indices.begin(); it != factory._loader->_instances[0]->B1[0].J_col_indices.end(); ++it) {
			dirichlet_indices.push_back(*it - 1);
		}
		dirichlet_values.insert(dirichlet_values.end(), factory._loader->_instances[0]->B1c[0].begin(), factory._loader->_instances[0]->B1c[0].end());
	}

	void fillNeighbours(std::vector<FETI4IMPIInt> &neighbours)
	{
		neighbours = factory._mesh->neighbours();
	}

protected:
	ECFConfiguration ecf;
	Factory factory;
	Step step;
};

}

int main(int argc, char** argv)
{
	// Always initialize MPI before call ESPRESO!!
	MPI_Init(&argc, &argv);

	espreso::APITestESPRESODataProvider provider(&argc, &argv);

	FETI4IMatrix K;
	FETI4IInt    matrixType = provider.matrixType();
	FETI4IInt    indexBase = 0;

	FETI4ICreateStiffnessMatrix(&K, matrixType, indexBase);

	for (size_t e = 0; e < provider.elements(); e++) {
		provider.addElementMatrix(K, e);
	}

	FETI4IInt iopts[FETI4I_INTEGER_OPTIONS_SIZE];
	FETI4IReal ropts[FETI4I_REAL_OPTIONS_SIZE];

	FETI4ISetDefaultIntegerOptions(iopts);
	FETI4ISetDefaultRealOptions(ropts);

	FETI4IInstance            instance;
	std::vector<FETI4IReal>   rhs;
	std::vector<FETI4IInt>    dirichlet_indices;
	std::vector<FETI4IReal>   dirichlet_values;
	std::vector<FETI4IInt>    l2g;
	std::vector<FETI4IMPIInt> neighbours;

	provider.computeRHS(rhs);
	provider.fillL2G(l2g);
	provider.fillDirichlet(dirichlet_indices, dirichlet_values);
	provider.fillNeighbours(neighbours);

	FETI4ICreateInstance(
			&instance,
			K, // Stiffness matrix
			rhs.size(), rhs.data(), // RHS
			l2g.data(), // local 2 global indices mapping
			neighbours.size(), neighbours.data(), // neighbours clusters
			dirichlet_indices.size(), dirichlet_indices.data(), dirichlet_values.data(), // Dirichlet boundary condition
			iopts, ropts); // FETI4I options

	std::vector<FETI4IReal> solution(rhs.size());

	FETI4ISolve(instance, solution.size(), solution.data());

	FETI4IDestroy(K);
	FETI4IDestroy(instance);

	MPI_Finalize();
}

//int main(int argc, char** argv)
//{
//	// Always initialize MPI before call ESPRESO!
//	MPI_Init(&argc, &argv);
//
/////////////////////////////////////////////// CREATE INSTANCE /////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	espreso::GlobalConfiguration configuration(&argc, &argv);
//
//	ESINFO(espreso::OVERVIEW) << "Run ESPRESO API test on " << espreso::environment->MPIsize << " process(es).";
//
//	// use ESPRESO factory to allow run test on all examples/assemblers
//	espreso::Factory factory(configuration);
//	factory.mesh->partitiate(1); // we need only one sub-domain per cluster
//	factory.instance->init(); // it perform full initialization; it is not effective but sufficient to test API
//
//	std::vector<espreso::Property> DOFs = factory.instance->physics().pointDOFs;
//	size_t DOFsSize = factory.instance->physics().f[0].size();
//
//	std::vector<FETI4IInt> l2g(DOFsSize);
//
//	// TODO: generalize l2g and dirichlet
//	for (size_t n = 0; n < factory.mesh->coordinates().clusterSize(); n++) {
//		for (size_t dof = 0; dof < DOFs.size(); dof++) {
//			l2g[DOFs.size() * n + dof] = DOFs.size() * factory.mesh->coordinates().globalIndex(n) + dof;
//		}
//	}
//
//	std::vector<FETI4IInt> dirichletIndices;
//	std::vector<FETI4IReal> dirichletValues;
//	for (size_t n = 0; n < factory.mesh->nodes().size(); n++) {
//		for (size_t dof = 0; dof < DOFs.size(); dof++) {
//			if (factory.mesh->nodes()[n]->hasProperty(DOFs[dof], 0)) {
//				dirichletIndices.push_back(DOFs.size() * n + dof);
//				dirichletValues.push_back(factory.mesh->nodes()[n]->getProperty(DOFs[dof], 0, 0, 0));
//			}
//		}
//	}
///////////////////////////////////////////////// USE ESPRESO API /////////////////////////////////////////////////////////
//
//	// At first create stiffness matrix
//	FETI4IMatrix K;
//	FETI4IInt indexBase = 0;
//	FETI4ICreateStiffnessMatrix(&K, (int)factory.instance->physics().mtype, indexBase);
//
//	std::vector<FETI4IReal> rhs(factory.instance->physics().f[0]);
//	std::vector<FETI4IMPIInt> neighbours(factory.mesh->neighbours());
//
//	// Compose the matrix from elements matrices
//	for (size_t e = 0; e < factory.mesh->elements().size(); e++) {
//		FETI4IInt dimension = static_cast<FETI4IInt>(factory.mesh->elements()[e]->type());
//		std::vector<FETI4IInt> nodes;
//		std::vector<FETI4IInt> dofs;
//		espreso::DenseMatrix Ke;
//		std::vector<double> fe;
//
//		for (size_t n = 0; n < factory.mesh->elements()[e]->nodes(); n++) {
//			nodes.push_back(factory.mesh->elements()[e]->node(n));
//		}
//
//		factory.instance->physics().assembleStiffnessMatrix(factory.mesh->elements()[e], Ke, fe, dofs);
//		FETI4IAddElement(K, dimension, nodes.size(), nodes.data(), dofs.size(), dofs.data(), Ke.values());
//	}
//
//	FETI4IInt iopts[FETI4I_INTEGER_OPTIONS_SIZE];
//	FETI4IReal ropts[FETI4I_REAL_OPTIONS_SIZE];
//
//	// FETI4ISetDefaultIntegerOptions(iopts);
//	// FETI4ISetDefaultRealOptions(ropts);
//
//	// Configure ESPRESO solver
//	iopts[FETI4I_FETI_METHOD] = static_cast<int>(configuration.linear_elasticity_2D.espreso.method);
//	iopts[FETI4I_CGSOLVER] = static_cast<int>(configuration.linear_elasticity_3D.espreso.solver);
//	iopts[FETI4I_SUBDOMAINS] = configuration.esdata.domains;
//	iopts[FETI4I_MAX_ITERATIONS] = configuration.linear_elasticity_3D.espreso.iterations;
//	iopts[FETI4I_PRECONDITIONER] = static_cast<int>(configuration.linear_elasticity_3D.espreso.preconditioner);
//	iopts[FETI4I_VERBOSE_LEVEL] = configuration.env.verbose_level;
//	iopts[FETI4I_TESTING_LEVEL] = configuration.env.testing_level;
//	iopts[FETI4I_MEASURE_LEVEL] = configuration.env.measure_level;
//	iopts[FETI4I_PRINT_MATRICES] = configuration.env.print_matrices;
//	ropts[FETI4I_PRECISION] = configuration.linear_elasticity_3D.espreso.epsilon;
//
//
//	// Create instance of a problem
//	FETI4IInstance instance;
//	FETI4ICreateInstance(
//			&instance,
//			K,
//			rhs.size(),
//			rhs.data(),
//			l2g.data(),
//			neighbours.size(),
//			neighbours.data(),
//			dirichletIndices.size(),
//			dirichletIndices.data(),
//			dirichletValues.data(),
//			iopts,
//			ropts);
//
//	// Prepare memory for save solution
//	std::vector<std::vector<FETI4IReal> > solution(1, std::vector<FETI4IReal>(rhs.size()));
//
//	// Solve the system
//	FETI4ISolve(instance, solution[0].size(), solution[0].data());
//
//	// Process solution
//
//	espreso::output::VTKXMLASCII vtk(configuration.output, factory.mesh, "results");
//	espreso::Step step;
//	std::vector<espreso::Solution*> sol;
//	sol.push_back(new espreso::Solution("solution", espreso::ElementType::NODES, factory.instance->physics().pointDOFs, solution));
//	vtk.storeSolution(step, sol);
//
//	// Remove data
//	FETI4IDestroy(K);
//	FETI4IDestroy(instance);
//
//	MPI_Finalize();
//}



