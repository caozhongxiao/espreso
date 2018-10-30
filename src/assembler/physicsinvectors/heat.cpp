
#include "heat.h"

#include "../instance.h"

#include "../../basis/containers/serializededata.h"
#include "../../basis/evaluator/evaluator.h"
#include "../../basis/matrices/denseMatrix.h"

#include "../../mesh/mesh.h"
#include "../../mesh/elements/element.h"
#include "../../mesh/store/elementstore.h"
#include "../../mesh/store/nodestore.h"

#include "../../config/ecf/environment.h"
#include "../../config/ecf/physics/heattransfer.h"

using namespace espreso;

Heat::Heat(Mesh &mesh, Instance &instance, Step &step, const HeatTransferConfiguration &configuration)
: PhysicsInVectors("HEAT", mesh, instance, step, configuration),
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

void Heat::processElement(eslocal domain, eslocal eindex, DenseMatrix &Ke, DenseMatrix &fe)
{
	auto epointer = _mesh.elements->epointers->datatarray()[eindex];

	const std::vector<DenseMatrix> &N = *(epointer->N);
	const std::vector<DenseMatrix> &dN = *(epointer->dN);
	const std::vector<double> &weighFactor = *(epointer->weighFactor);

	DenseMatrix Ce(2, 2), coordinates(epointer->nodes, 2), J(2, 2), invJ(2, 2), K(epointer->nodes, 4), dND;
	double detJ;
	DenseMatrix gpK(1, 4);

	for (int n = 0; n < epointer->nodes; n++) {
		coordinates(n, 0) = (_coordinates->begin() + n)->at(0);
		coordinates(n, 1) = (_coordinates->begin() + n)->at(1);
		K(n, 0) = (_K->begin() + n)->at(0);
		K(n, 1) = (_K->begin() + n)->at(1);
		K(n, 2) = (_K->begin() + n)->at(2);
		K(n, 3) = (_K->begin() + n)->at(3);
	}

	eslocal Ksize = epointer->nodes;

	Ke.resize(Ksize, Ksize);
	fe.resize(Ksize, 1);
	Ke = 0;
	fe = 0;

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

		for (eslocal i = 0; i < Ksize; i++) {
			fe(i, 0) += detJ * weighFactor[gp] * N[gp](0, i);
		}
	}

	std::cout << Ke << fe;
}
