
#include "laplacesteklovpoincare3d.h"
#include "../../config/ecf/physics/heattransfer.h"

#include "../../basis/logging/logging.h"
#include "../../basis/matrices/sparseVVPMatrix.h"

#include "../instance.h"
#include "../step.h"

#include "../../old/mesh/structures/mesh.h"
#include "../../old/mesh/structures/coordinates.h"
#include "../../old/mesh/structures/elementtypes.h"

#include "../../solver/generic/SparseMatrix.h"

#ifdef BEM4I
#include "esbem.h"
#endif

#include "../../basis/utilities/utils.h"

using namespace espreso;

LaplaceSteklovPoincare3D::LaplaceSteklovPoincare3D(Mesh *mesh, Instance *instance, const HeatTransferConfiguration &configuration, const ResultsSelectionConfiguration &propertiesConfiguration)
: Physics("LAPLACE STEKLOV POINCARE 3D", mesh, instance, &configuration), HeatTransfer3D(mesh, instance, configuration, propertiesConfiguration)
{
#ifndef BEM4I
	ESINFO(GLOBAL_ERROR) << "BEM4I is not linked!. Copy BEM4I library to tools/bem4i and re-configure ESPRESO.";
#endif

	for (size_t loadStep = 0; loadStep < configuration.load_steps; loadStep++) {
		if (configuration.load_steps_settings.at(loadStep + 1).convection.size() || configuration.load_steps_settings.at(loadStep + 1).diffuse_radiation.size()) {
			ESINFO(GLOBAL_ERROR) << "BEM discretization not supports CONVECTION or RADIATION boundary condition.";
		}
	}
}

void LaplaceSteklovPoincare3D::prepareHybridTotalFETIWithKernels()
{
	// extraction of boundary nodes compute all faces on domains. Hence, no face computation is needed
	prepare();
}

void LaplaceSteklovPoincare3D::preprocessData(const Step &step)
{
	// TODO: MESH
//	if (offset == (size_t)-1) {
//		offset = _instance->solutions.size();
//		_instance->solutions.resize(offset + SolutionIndex::SIZE, NULL);
//
//		_instance->solutions[offset + SolutionIndex::TEMPERATURE] = new Solution(*_mesh, "temperature", ElementType::NODES, 1);
//		computeInitialTemperature(step, _instance->solutions[offset + SolutionIndex::TEMPERATURE]->data);
//	}
//
//	computeInitialTemperature(step, _instance->primalSolution);
//	temperature = _mesh->nodes->appendData({ "TEMPERATURE" }, _instance->primalSolution);
//
//	if (BEMOffset == (size_t)-1) {
//		BEMOffset = _instance->solutions.size();
//		_instance->solutions.resize(BEMOffset + BEMSolutionIndex::SIZE, NULL);
//		_instance->solutions[BEMOffset + BEMSolutionIndex::BOUNDARY] = new Solution(*_mesh, "temperatureOnBoundary", ElementType::NODES, 1, _instance->primalSolution);
//	}
}

void LaplaceSteklovPoincare3D::updateMatrix(const Step &step, Matrices matrices, size_t domain)
{
	if (matrices & Matrices::f) {
		_instance->f[domain].clear();
		_instance->f[domain].resize(_instance->domainDOFCount[domain]);
	}

	std::vector<eslocal> elements;
	std::vector<double> coordinates;
	// TODO: MESH
//	boundaryTriangularization(elements, coordinates, domain);
//
//	_instance->K[domain].rows = _boundaryIndices[domain].size();
//	_instance->K[domain].cols = _boundaryIndices[domain].size();
//	_instance->K[domain].nnz  = _boundaryIndices[domain].size() * _boundaryIndices[domain].size();
//	_instance->K[domain].type = 'G';
//	_instance->K[domain].dense_values.resize(_instance->K[domain].nnz);
//	_instance->K[domain].mtype = MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;

#ifdef BEM4I
//	bem4i::getLaplaceSteklovPoincare(
//			_instance->K[domain].dense_values.data(),
//			(eslocal)_boundaryIndices[domain].size(),
//			coordinates.data(),
//			(eslocal)(elements.size() / 3),
//			elements.data(),
//			3, 3, 0);
#endif

	_instance->K[domain].ConvertDenseToCSR(1);
	SparseVVPMatrix<eslocal> K, M; // never filled
	assembleBoundaryConditions(K, M, step, Matrices::f, domain);
}

void LaplaceSteklovPoincare3D::updateMatrix(const Step &step, Matrices matrices, eslocal eindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe)
{
	ESINFO(ERROR) << "BEM discretization not supports assembling of stiffness matrix K for one element.";
}

void LaplaceSteklovPoincare3D::processSolution(const Step &step)
{
	// TODO: get solution for all nodes from BEM library
//	#pragma omp parallel for
//	for (size_t p = 0; p < _mesh->parts(); p++) {
//		std::fill(_instance->solutions[offset + SolutionIndex::TEMPERATURE]->data[p].begin(), _instance->solutions[offset + SolutionIndex::TEMPERATURE]->data[p].end(), 0);
//		for (size_t i = 0; i < _boundaryIndices[p].size(); i++) {
//			_instance->solutions[offset + SolutionIndex::TEMPERATURE]->data[p][_mesh->coordinates().localIndex(_boundaryIndices[p][i], p)] = _instance->primalSolution[p][i];
//		}
//	}
}



