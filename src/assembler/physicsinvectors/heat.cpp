
#include "heat.h"

#include "../instance.h"

#include "../../basis/containers/serializededata.h"
#include "../../basis/evaluator/evaluator.h"
#include "../../basis/matrices/denseMatrix.h"
#include "../../basis/utilities/utils.h"

#include "../../mesh/mesh.h"
#include "../../mesh/elements/element.h"
#include "../../mesh/store/elementstore.h"
#include "../../mesh/store/nodestore.h"
#include "../../mesh/store/boundaryregionstore.h"

#include "../../config/ecf/environment.h"
#include "../../config/ecf/physics/heattransfer.h"

#include "../../solver/generic/SparseMatrix.h"

using namespace espreso;

Heat::Heat(Mesh &mesh, Instance &instance, Step &step, const HeatTransferConfiguration &configuration)
: PhysicsInVectors("HEAT", mesh, instance, step, configuration),
  _configuration(configuration),
  _coordinates(NULL), _K(NULL),
  _area(NULL),
  _temperature(NULL)
{
	_coordinates = new serializededata<eslocal, double>(2, _nDistribution);
//	_temperature = new serializededata<eslocal, double>(1, _nDistribution);
	_K = new serializededata<eslocal, double>(4, _nDistribution);

	_area = new serializededata<eslocal, double>(1, _nDistribution);

	_temperature = _mesh.nodes->appendData(1, { }, _instance.primalSolution);
}

void Heat::initData()
{
	size_t threads = environment->OMP_NUM_THREADS;

	const MaterialConfiguration* material = _mesh.materials.front();
	tarray<double> kxx(_mesh.nodes->distribution, 1);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		material->thermal_conductivity.values.get(0, 0).evaluator->evaluate(kxx.size(t), NULL, NULL, 0, kxx.begin(t));

		auto c = _coordinates->begin(t);
		auto K = _K->begin(t);
		for (auto n = _mesh.elements->procNodes->datatarray().begin(t); n != _mesh.elements->procNodes->datatarray().end(t); ++n, ++c, ++K) {
			c->at(0) = _mesh.nodes->coordinates->datatarray()[*n].x;
			c->at(1) = _mesh.nodes->coordinates->datatarray()[*n].y;
			K->at(0) = kxx[*n];
			K->at(1) = kxx[*n];
		}
	}

//	std::cout << "Coordinates: " << *_coordinates << "\n";
//	std::cout << "Conductivity: " << *_K << "\n";
}

void Heat::updateData()
{
	// update data that could be changed
}

static inline double determinant2x2(double *values)
{
	double det = values[0] * values[3] - values[1] * values[2];
	if (det < 0) {
		printf("negative determinant\n");
		exit(0);
	}
	return det;
}

static inline void inverse2x2(const double *m, double *inv, double det)
{
	double detJx = 1 / det;
	inv[0] =   detJx * m[3];
	inv[1] = - detJx * m[1];
	inv[2] = - detJx * m[2];
	inv[3] =   detJx * m[0];
}

eslocal Heat::processElement(eslocal eindex, eslocal nindex, DenseMatrix &Ke, DenseMatrix &fe)
{
	auto epointer = _mesh.elements->epointers->datatarray()[eindex];

	const std::vector<DenseMatrix> &N = *(epointer->N);
	const std::vector<DenseMatrix> &dN = *(epointer->dN);
	const std::vector<double> &weighFactor = *(epointer->weighFactor);

	DenseMatrix Ce(2, 2), coordinates(epointer->nodes, 2), J(2, 2), invJ(2, 2), K(epointer->nodes, 4), dND;
	double detJ;
	DenseMatrix gpK(1, 4);

	for (int n = 0; n < epointer->nodes; n++) {
		coordinates(n, 0) = (_coordinates->begin() + nindex + n)->at(0);
		coordinates(n, 1) = (_coordinates->begin() + nindex + n)->at(1);
		K(n, 0) = (_K->begin() + nindex + n)->at(0);
		K(n, 1) = (_K->begin() + nindex + n)->at(1);
		K(n, 2) = (_K->begin() + nindex + n)->at(2);
		K(n, 3) = (_K->begin() + nindex + n)->at(3);
	}

	eslocal Ksize = epointer->nodes;

	Ke.resize(Ksize, Ksize);
	fe.resize(Ksize, 1);
	Ke = 0;
	fe = 0;

//	if (eindex == 0) {
//		std::cout << "C: " << coordinates;
//		std::cout << "K: " << K;
//	}

	for (size_t gp = 0; gp < N.size(); gp++) {
		J.multiply(dN[gp], coordinates);
		detJ = determinant2x2(J.values());
		inverse2x2(J.values(), invJ.values(), detJ);

		gpK.multiply(N[gp], K);

		Ce(0, 0) = gpK(0, 0);
		Ce(1, 1) = gpK(0, 1);
		Ce(0, 1) = gpK(0, 2);
		Ce(1, 0) = gpK(0, 3);

		dND.multiply(invJ, dN[gp]);

		Ke.multiply(dND, Ce * dND, detJ * weighFactor[gp], 1, true);
	}

//	if (eindex == 0) {
//		std::cout << " xx " << eindex << " xx \n" << Ke << fe;
//	}
	return epointer->nodes;
}

void Heat::setDirichlet()
{
//	std::cout << _instance.K.front();
//	std::cout << _instance.f.front();

	eslocal halosize = _instance.K.front().haloRows;

	auto dirichlet = _configuration.load_steps_settings.at(1).temperature;
	for (auto it = dirichlet.regions.begin(); it != dirichlet.regions.end(); ++it) {
		BoundaryRegionStore *region = _mesh.bregion(it->first);
		ECFExpression &expression = it->second;

		auto &ROW = _instance.K.front().CSR_I_row_indices;
		auto &COL = _instance.K.front().CSR_J_col_indices;
		auto &VAL = _instance.K.front().CSR_V_values;
		auto &RHS = _instance.f.front();

		for (auto r = region->uniqueNodes->datatarray().begin(); r != region->uniqueNodes->datatarray().end(); ++r) {
//			std::cout << "r: " << *r << "\n";
			eslocal row = _globalIndices[*r] - _mesh.nodes->uniqueOffset + halosize;
			for (eslocal i = ROW[row]; i < ROW[row + 1]; i++) {
				RHS[row] = expression.evaluator->evaluate(_mesh.nodes->coordinates->datatarray()[*r], 0, 0);
				if (COL[i - 1] - 1 == _globalIndices[*r]) {
					VAL[i - 1] = 1;
				} else {
					VAL[i - 1] = 0;
					eslocal col = std::lower_bound(_globalIndices.begin(), _globalIndices.end(), COL[i - 1] - 1) - _globalIndices.begin();
//					printf("%d col: %d -> %d\n", environment->MPIrank, COL[i - 1], col);
					for (eslocal c = ROW[col]; c < ROW[col + 1]; c++) {
						if (COL[c - 1] - 1 == _globalIndices[*r]) {
//							RHS[col] -= VAL[c - 1] * RHS[row];
							VAL[c - 1] = 0;
						}
					}
				}
			}
		}
	}

//	std::cout << _instance.K.front();
//	std::cout << _instance.f.front();
}
