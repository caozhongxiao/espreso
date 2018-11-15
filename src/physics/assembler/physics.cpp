
#include "physics.h"

#include "../../basis/containers/serializededata.h"
#include "../../basis/utilities/utils.h"
#include "../../basis/utilities/parser.h"
#include "../../basis/utilities/communication.h"
#include "../instance.h"
#include "../step.h"

#include "../constraints/constraints.h"

#include "../../mesh/mesh.h"
#include "../../mesh/elements/element.h"
#include "../../mesh/store/nodestore.h"
#include "../../mesh/store/elementstore.h"
#include "../../mesh/store/boundaryregionstore.h"
#include "../../mesh/store/elementsregionstore.h"

#include "../../solver/generic/SparseMatrix.h"
#include "../../basis/matrices/sparseVVPMatrix.h"
#include "../../basis/matrices/denseMatrix.h"
#include "../../basis/matrices/sparseCSRMatrix.h"
#include "../../config/ecf/environment.h"
#include "../../config/ecf/solver/feti.h"
#include "../../config/ecf/physics/physics.h"

#include <algorithm>
#include <numeric>

#ifdef BEM4I
#include "esbem.h"
#endif

using namespace espreso;

Physics::Physics()
: _DOF(NULL), _name(""), _mesh(NULL), _instance(NULL), _step(NULL), _constraints(NULL), _configuration(NULL), _DOFs(0), _invalidElements(0)
{

}

Physics::Physics(const std::string &name, Mesh *mesh, Instance *instance, Step *step, const PhysicsConfiguration *configuration, int DOFs)
: _DOF(NULL), _name(name), _mesh(mesh), _instance(instance), _step(step), _constraints(NULL), _configuration(configuration), _DOFs(DOFs), _invalidElements(0) // initialized in a particular physics
{
	std::vector<int> BEMRegions(_mesh->elements->regionMaskSize);
	for (auto it = configuration->discretization.begin(); it != configuration->discretization.end(); ++it) {
		if (it->second == DISCRETIZATION::BEM) {
			for (size_t r = 0; r < _mesh->elementsRegions.size(); r++) {
				if (StringCompare::caseInsensitiveEq(it->first, _mesh->elementsRegions[r]->name)) {
					BEMRegions[r / (8 * sizeof(int))] |= 1 << (r % (8 * sizeof(int)));
				}
			}
		}
	}

	_hasBEM = false;
	_BEMDomain.resize(_mesh->elements->ndomains);
	for (eslocal d = 0; d < _mesh->elements->ndomains; d++) {
		for (int i = 0; i < _mesh->elements->regionMaskSize; i++) {
			if (_mesh->elements->regions->datatarray()[_mesh->elements->elementsDistribution[d] * _mesh->elements->regionMaskSize + i] & BEMRegions[i]) {
				_BEMDomain[d] = 1;
				_hasBEM = true;
			}
		}
	}

	_BEMData.resize(mesh->elements->ndomains, NULL);
}

Physics::~Physics()
{
	if (_constraints != NULL) {
		delete _constraints;
	}
#ifdef BEM4I
	for (size_t i = 0; i < _BEMData.size(); i++) {
		if (_BEMData[i] != NULL) {
			bem4i::deleteBem4iData(_BEMData[i]);
		}
	}
#endif
}


void Physics::initLocalNodeUniformDOFs(std::vector<eslocal> &offsets, eslocal multiplier)
{
	size_t threads = environment->OMP_NUM_THREADS;

	// nID, domain
	std::vector<std::vector<std::pair<eslocal, eslocal> > > ntodomains(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (eslocal d = _mesh->elements->domainDistribution[t]; d != _mesh->elements->domainDistribution[t + 1]; ++d) {
			std::vector<eslocal> dnodes(
					(_mesh->elements->procNodes->begin() + _mesh->elements->elementsDistribution[d])->begin(),
					(_mesh->elements->procNodes->begin() + _mesh->elements->elementsDistribution[d + 1])->begin());

			Esutils::sortAndRemoveDuplicity(dnodes);
			for (size_t i = 0; i < dnodes.size(); i++) {
				ntodomains[t].push_back(std::pair<eslocal, eslocal>(dnodes[i], d));
			}
		}
	}

	for (size_t t = 1; t < threads; t++) {
		ntodomains[0].insert(ntodomains[0].end(), ntodomains[t].begin(), ntodomains[t].end());
	}

	std::sort(ntodomains[0].begin(), ntodomains[0].end());

	std::vector<std::vector<eslocal> > DOFs(threads);
	std::vector<size_t> ndistribution = tarray<size_t>::distribute(threads, ntodomains[0].size());

	for (size_t t = 1; t < threads; t++) {
		while (
				ndistribution[t] < ntodomains[0].size() &&
				ntodomains[0][ndistribution[t]].first == ntodomains[0][ndistribution[t] - 1].first) {

			++ndistribution[t];
		}
	}

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<eslocal> tDOFs(_mesh->elements->ndomains);

		for (size_t n = ndistribution[t]; n < ndistribution[t + 1]; ++n) {
			tDOFs[ntodomains[0][n].second] += multiplier;
		}

		DOFs[t].swap(tDOFs);
	}

	Esutils::sizesToOffsets(DOFs, offsets);

	std::vector<std::vector<std::vector<eslocal> > > sBuffer(threads);
	std::vector<std::vector<eslocal> > rBuffer(_mesh->neighboursWithMe.size());


	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto nranks = _mesh->nodes->ranks->begin() + ntodomains[0][ndistribution[t]].first;

		std::vector<std::vector<eslocal> > tBuffer(_mesh->neighboursWithMe.size());

		size_t n = ndistribution[t];
		while (n < ndistribution[t + 1]) {
			size_t begin = n++;
			while (n < ntodomains[0].size() && ntodomains[0][n].first == ntodomains[0][n - 1].first) {
				++n;
			}

			eslocal noffset = 0;
			for (auto r = nranks->begin(); r != nranks->end(); ++r) {
				while (_mesh->neighboursWithMe[noffset] < *r) {
					++noffset;
				}

				tBuffer[noffset].push_back(n - begin);
				for (size_t i = begin; i < n; i++) {
					tBuffer[noffset].push_back(_mesh->elements->firstDomain + ntodomains[0][i].second);
					for (size_t dof = 0; dof < multiplier; ++dof) {
						tBuffer[noffset].push_back(DOFs[t][ntodomains[0][i].second] + dof);
					}
				}
			}
			++nranks;

			for (size_t i = begin; i < n; i++) {
				DOFs[t][ntodomains[0][i].second] += multiplier;
			}
		}
		sBuffer[t].swap(tBuffer);
	}

	for (size_t t = 1; t < threads; t++) {
		for (size_t n = 0; n < sBuffer[t].size(); n++) {
			sBuffer[0][n].insert(sBuffer[0][n].end(), sBuffer[t][n].begin(), sBuffer[t][n].end());
		}
	}

	if (!Communication::exchangeUnknownSize(sBuffer[0], rBuffer, _mesh->neighboursWithMe)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange uniform DOFs.";
	}

	// domain0, DOF0, DOF1, ..., DOFn, domain1, DOF0, DOF1, ..., DOFn, ...; domain0, ...
	std::vector<eslocal> DOFDistribution(1), DOFData;

	// TODO: make it parallel
	// parallelization is possible if node order will be kept as: boundary first!
	// now we prefer generality
	auto nranks = _mesh->nodes->ranks->begin();
	std::vector<eslocal> roffset(rBuffer.size());
	for (eslocal n = 0; n < _mesh->nodes->size; ++n, ++nranks) {
		eslocal noffset = 0;
		for (auto r = nranks->begin(); r != nranks->end(); ++r) {
			while (_mesh->neighboursWithMe[noffset] < *r) {
				++noffset;
			}

			eslocal domains = rBuffer[noffset][roffset[noffset]++];
			for (eslocal d = 0; d < domains; ++d) {
				DOFData.push_back(rBuffer[noffset][roffset[noffset]++]);
				for (size_t dof = 0; dof < multiplier; ++dof) {
					DOFData.push_back(rBuffer[noffset][roffset[noffset]++]);
				}
			}
		}
		DOFDistribution.push_back(DOFData.size());
	}

	std::vector<size_t> distribution = _mesh->nodes->distribution, datadistribution(threads + 1);
	for (size_t t = 1; t < threads; t++) {
		++distribution[t];
		datadistribution[t] = DOFDistribution[distribution[t]];
	}
	datadistribution[threads] = DOFDistribution[distribution[threads]];
	++distribution[threads];


	_DOF = new serializededata<eslocal, eslocal>(
			tarray<eslocal>(distribution, DOFDistribution),
			tarray<eslocal>(datadistribution, DOFData));

//	Communication::serialize([&] () {
//		std::cout << " -- " << environment->MPIrank << " -- \n";
//		std::cout << *_DOF << "\n";
//	});
}

void Physics::initGlobalNodeUniformDOFs(eslocal &offset, eslocal multiplier)
{

}

struct IJ {
	eslocal row, column;
};

bool operator==(const IJ &left, const IJ &right)
{
	return left.row == right.row && left.column == right.column;
}

bool operator!=(const IJ &left, const IJ &right)
{
	return !(left == right);
}

bool operator<(const IJ &left, const IJ &right)
{
	return left.row == right.row ? left.column < right.column : left.row < right.row;
}

void Physics::buildLocalNodeUniformCSRPattern(eslocal multiplier)
{
	size_t threads = environment->OMP_NUM_THREADS;

	_KPermutation.resize(_mesh->elements->ndomains);
	_RHSPermutation.resize(_mesh->elements->ndomains);

	_instance->K.resize(_mesh->elements->ndomains);
	_instance->M.resize(_mesh->elements->ndomains);
	_instance->f.resize(_mesh->elements->ndomains);
	_instance->R.resize(_mesh->elements->ndomains);
	_instance->primalSolution.resize(_mesh->elements->ndomains);
	for (eslocal d = 0; d < _mesh->elements->ndomains; d++) {
		_instance->K[d].rows = _instance->domainDOFCount[d];
		_instance->K[d].cols = _instance->domainDOFCount[d];
		_instance->f[d].resize(_instance->domainDOFCount[d]);
		_instance->R[d].resize(_instance->domainDOFCount[d]);
	}

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (eslocal d = _mesh->elements->domainDistribution[t]; d != _mesh->elements->domainDistribution[t + 1]; ++d) {
			std::vector<eslocal> permK, permRHS;
			std::vector<IJ> KPattern;
			std::vector<eslocal> RHSPattern, ROW, COL;

			auto fullInsert = [&] (serializededata<eslocal, eslocal>::const_iterator &enodes) {
				size_t size = RHSPattern.size();
				for (auto n = enodes->begin(); n != enodes->end(); ++n) {
					auto DOFs = (_DOF->begin() + *n)->begin();
					while (*DOFs != d + _mesh->elements->firstDomain) {
						DOFs += 1 + multiplier;
					}
					RHSPattern.insert(RHSPattern.end(), DOFs + 1, DOFs + 1 + multiplier);
				}
				for (auto row = RHSPattern.begin() + size; row != RHSPattern.end(); ++row) {
					for (auto col = RHSPattern.begin() + size; col != RHSPattern.end(); ++col) {
						KPattern.push_back({*row, *col});
					}
				}
			};

			auto upperInsert = [&] (serializededata<eslocal, eslocal>::const_iterator &enodes) {
				size_t size = RHSPattern.size();
				for (auto n = enodes->begin(); n != enodes->end(); ++n) {
					auto DOFs = (_DOF->begin() + *n)->begin();
					while (*DOFs != d + _mesh->elements->firstDomain) {
						DOFs += 1 + multiplier;
					}
					RHSPattern.insert(RHSPattern.end(), DOFs + 1, DOFs + 1 + multiplier);
				}
				for (auto row = RHSPattern.begin() + size, colbegin = RHSPattern.begin() + size; row != RHSPattern.end(); ++row, ++colbegin) {
					for (auto col = colbegin; col != RHSPattern.end(); ++col) {
						KPattern.push_back(*row <= *col ? IJ{*row, *col} : IJ{*col, *row});
					}
				}
			};

			auto ebegin = _mesh->elements->procNodes->cbegin() + _mesh->elements->elementsDistribution[d];
			auto eend = _mesh->elements->procNodes->cbegin() + _mesh->elements->elementsDistribution[d + 1];


			for (auto e = ebegin; e != eend; ++e) {
				switch (getMatrixType(d)) {
				case MatrixType::REAL_UNSYMMETRIC: fullInsert(e); break;
				default: upperInsert(e);
				}
			}

			for (size_t r = 0; r < _mesh->boundaryRegions.size(); r++) {
				if (
						_mesh->boundaryRegions[r]->dimension &&
						_mesh->boundaryRegions[r]->eintervalsDistribution[d] < _mesh->boundaryRegions[r]->eintervalsDistribution[d + 1]) {

					eslocal begin = _mesh->boundaryRegions[r]->eintervals[_mesh->boundaryRegions[r]->eintervalsDistribution[d]].begin;
					eslocal end = _mesh->boundaryRegions[r]->eintervals[_mesh->boundaryRegions[r]->eintervalsDistribution[d + 1] - 1].end;
					auto enodes = _mesh->boundaryRegions[r]->procNodes->cbegin() + begin;

					for (eslocal i = begin; i < end; ++i, ++enodes) {
						switch (getMatrixType(d)) {
						case MatrixType::REAL_UNSYMMETRIC: fullInsert(enodes); break;
						default: upperInsert(enodes);
						}
					}
				}
			}

			std::vector<eslocal> pK(KPattern.size());
			std::iota(pK.begin(), pK.end(), 0);
			std::sort(pK.begin(), pK.end(), [&] (eslocal i, eslocal j) {
				return KPattern[i] < KPattern[j];
			});

			std::vector<eslocal> pRHS(RHSPattern.size());
			std::iota(pRHS.begin(), pRHS.end(), 0);
			std::sort(pRHS.begin(), pRHS.end(), [&] (eslocal i, eslocal j) {
				return RHSPattern[i] < RHSPattern[j];
			});

			ROW.reserve(_instance->domainDOFCount[d] + 1);
			COL.reserve(KPattern.size());
			ROW.push_back(1);
			COL.push_back(KPattern[pK.front()].column + 1);
			permK.resize(KPattern.size());
			permK[pK.front()] = 0;
			for (size_t i = 1, nonzeros = 0; i < pK.size(); i++) {
				if (KPattern[pK[i]] != KPattern[pK[i - 1]]) {
					++nonzeros;
					COL.push_back(KPattern[pK[i]].column + 1);
					if (KPattern[pK[i - 1]].row != KPattern[pK[i]].row) {
						ROW.push_back(nonzeros + 1);
					}
				}
				permK[pK[i]] = nonzeros;
			}
			permRHS.resize(RHSPattern.size());
			permRHS[pRHS.front()] = 0;
			for (size_t i = 1, nonzeros = 0; i < pRHS.size(); i++) {
				if (RHSPattern[pRHS[i]] != RHSPattern[pRHS[i - 1]]) {
					++nonzeros;
				}
				permRHS[pRHS[i]] = nonzeros;
			}

			_KPermutation[d].swap(permK);
			_RHSPermutation[d].swap(permRHS);

			ROW.push_back(COL.size() + 1);

			_instance->K[d].nnz = COL.size();
			_instance->K[d].mtype = getMatrixType(d);
			switch (_instance->K[d].mtype) {
			case MatrixType::REAL_UNSYMMETRIC: _instance->K[d].type = 'G'; break;
			default: _instance->K[d].type = 'S';
			}
			_instance->K[d].CSR_V_values.resize(COL.size());
			_instance->K[d].CSR_I_row_indices.swap(ROW);
			_instance->K[d].CSR_J_col_indices.swap(COL);
			_instance->M[d] = _instance->K[d];
			_instance->M[d].mtype = MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;
			_instance->M[d].type = 'S';
		}
	}

//	std::cout << _KPermutation.front();
//	std::cout << _RHSPermutation.front();
}

void Physics::buildGlobalNodeUniformCSRPattern(eslocal multiplier)
{

}

void Physics::printInvalidElement(eslocal eindex) const
{
	if (_invalidElements++ == 0) {
		auto nodes = _mesh->elements->procNodes->begin() + eindex;

		std::ofstream os("invalidElement.vtk");
		os << "# vtk DataFile Version 2.0\n";
		os << "INVALID ELEMENT\n";
		os << "ASCII\n";
		os << "DATASET UNSTRUCTURED_GRID\n";
		os << "\n";
		os << "POINTS " << nodes->size() << " float\n";
		for (auto n = nodes->begin(); n != nodes->end(); ++n) {
			os << _mesh->nodes->coordinates->datatarray()[*n].x << " " << _mesh->nodes->coordinates->datatarray()[*n].y << " " << _mesh->nodes->coordinates->datatarray()[*n].z << "\n";
		}
		os << "\n";
		os << "CELLS 1 " << nodes->size() + 1 << "\n";
		os << nodes->size();
		for (auto n = nodes->begin(); n != nodes->end(); ++n) {
			os << " " << n - nodes->begin();
		}
		os << "\n";
		os << "CELL_TYPES 1\n";
		switch (_mesh->elements->epointers->datatarray()[eindex]->code) {
		case Element::CODE::SQUARE4:
			os << "9\n";
			break;
		case Element::CODE::SQUARE8:
			os << "23\n";
			break;
		case Element::CODE::TRIANGLE3:
			os << "5\n";
			break;
		case Element::CODE::TRIANGLE6:
			os << "22\n";
			break;
		case Element::CODE::TETRA4:
			os << "10\n";
			break;
		case Element::CODE::TETRA10:
			os << "24\n";
			break;
		case Element::CODE::PYRAMID5:
			os << "14\n";
			break;
		case Element::CODE::PYRAMID13:
			os << "27\n";
			break;
		case Element::CODE::PRISMA6:
			os << "13\n";
			break;
		case Element::CODE::PRISMA15:
			os << "26\n";
			break;
		case Element::CODE::HEXA8:
			os << "12\n";
			break;
		case Element::CODE::HEXA20:
			os << "25\n";
			break;
		default:
			break;
		}
	}
}

void Physics::updateMatrix(Matrices matrix)
{
	#pragma omp parallel for
	for  (size_t d = 0; d < _instance->domains; d++) {

		updateMatrix(matrix, d);

//		switch (_instance->K[d].mtype) {
//		case MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE:
//		case MatrixType::REAL_SYMMETRIC_INDEFINITE:
//			_instance->K[d].RemoveLower();
//			_instance->M[d].RemoveLower();
//			break;
//		case MatrixType::REAL_UNSYMMETRIC:
//			break;
//		}

		ESINFO(PROGRESS3) << Info::plain() << ".";
	}
	ESINFO(PROGRESS3);

	int my = _invalidElements, all = 0;
	MPI_Reduce(&my, &all, 1, MPI_INT, MPI_SUM, 0, environment->MPICommunicator);
	if (all) {
		ESINFO(ALWAYS_ON_ROOT) << Info::TextColor::YELLOW << "ESPRESO internal error: " << all << " invalid (dej <= 0) elements founded.";
	}
}

void Physics::updateMatrix(Matrices matrices, size_t domain)
{
	DenseMatrix Ke, Me, Re, fe;
	eslocal KIndex = 0, RHSIndex = 0;

	if (matrices & Matrices::K) {
		std::fill(_instance->K[domain].CSR_V_values.begin(), _instance->K[domain].CSR_V_values.end(), 0);
	}
	if (matrices & Matrices::M) {
		std::fill(_instance->M[domain].CSR_V_values.begin(), _instance->M[domain].CSR_V_values.end(), 0);
	}
	if (matrices & Matrices::f) {
		std::fill(_instance->f[domain].begin(), _instance->f[domain].end(), 0);
	}
	if (matrices & Matrices::R) {
		std::fill(_instance->R[domain].begin(), _instance->R[domain].end(), 0);
	}

	auto fullInsert = [&] (serializededata<eslocal, eslocal>::const_iterator &enodes, const double &Kreduction) {
		for (auto r = 0; r < enodes->size(); ++r, ++RHSIndex) {
			if ((matrices & Matrices::f) && fe.rows()) {
				_instance->f[domain][_RHSPermutation[domain][RHSIndex]] += _step->internalForceReduction * fe(r, 0);
			}
			if ((matrices & Matrices::R) && Re.rows()) {
				_instance->R[domain][_RHSPermutation[domain][RHSIndex]] += Re(r, 0);
			}

			for (auto c = 0; c < enodes->size(); ++c, ++KIndex) {
				if ((matrices & Matrices::K) && Ke.rows()) {
					_instance->K[domain].CSR_V_values[_KPermutation[domain][KIndex]] += Kreduction * Ke(r, c);
				}
				if ((matrices & Matrices::M) && Me.rows()) {
					_instance->M[domain].CSR_V_values[_KPermutation[domain][KIndex]] += Me(r, c);
				}
			}
		}
	};

	auto upperInsert = [&] (serializededata<eslocal, eslocal>::const_iterator &enodes, const double &Kreduction) {
		for (auto r = 0; r < enodes->size(); ++r, ++RHSIndex) {
			if ((matrices & Matrices::f) && fe.rows()) {
				_instance->f[domain][_RHSPermutation[domain][RHSIndex]] += _step->internalForceReduction * fe(r, 0);
			}
			if ((matrices & Matrices::R) && Re.rows()) {
				_instance->R[domain][_RHSPermutation[domain][RHSIndex]] += Re(r, 0);
			}

			for (auto c = r; c < enodes->size(); ++c, ++KIndex) {
				if ((matrices & Matrices::K) && Ke.rows()) {
					_instance->K[domain].CSR_V_values[_KPermutation[domain][KIndex]] += Kreduction * Ke(r, c);
				}
				if ((matrices & Matrices::M) && Me.rows()) {
					_instance->M[domain].CSR_V_values[_KPermutation[domain][KIndex]] += Me(r, c);
				}
			}
		}
	};

	auto boundary = [&] (Matrices restriction) {
		for (size_t r = 0; r < _mesh->boundaryRegions.size(); r++) {
			if (_mesh->boundaryRegions[r]->dimension == 2) {
				if (_mesh->boundaryRegions[r]->eintervalsDistribution[domain] < _mesh->boundaryRegions[r]->eintervalsDistribution[domain + 1]) {
					eslocal begin = _mesh->boundaryRegions[r]->eintervals[_mesh->boundaryRegions[r]->eintervalsDistribution[domain]].begin;
					eslocal end = _mesh->boundaryRegions[r]->eintervals[_mesh->boundaryRegions[r]->eintervalsDistribution[domain + 1] - 1].end;
					auto nodes = _mesh->boundaryRegions[r]->procNodes->cbegin() + begin;
					for (eslocal i = begin; i < end; ++i, ++nodes) {
						processFace(domain, _mesh->boundaryRegions[r], restriction, i, Ke, Me, Re, fe);
						switch (getMatrixType(domain)) {
						case MatrixType::REAL_UNSYMMETRIC: fullInsert(nodes, _step->internalForceReduction); break;
						default: upperInsert(nodes, _step->internalForceReduction);
						}
					}
				}
			}
			if (_mesh->boundaryRegions[r]->dimension == 1) {
				if (_mesh->boundaryRegions[r]->eintervalsDistribution[domain] < _mesh->boundaryRegions[r]->eintervalsDistribution[domain + 1]) {
					eslocal begin = _mesh->boundaryRegions[r]->eintervals[_mesh->boundaryRegions[r]->eintervalsDistribution[domain]].begin;
					eslocal end = _mesh->boundaryRegions[r]->eintervals[_mesh->boundaryRegions[r]->eintervalsDistribution[domain + 1] - 1].end;
					auto nodes = _mesh->boundaryRegions[r]->procNodes->cbegin() + begin;
					for (eslocal i = begin; i < end; ++i, ++nodes) {
						processEdge(domain, _mesh->boundaryRegions[r], restriction, i, Ke, Me, Re, fe);
						switch (getMatrixType(domain)) {
						case MatrixType::REAL_UNSYMMETRIC: fullInsert(nodes, _step->internalForceReduction); break;
						default: upperInsert(nodes, _step->internalForceReduction);
						}
					}
				}
			}
		}
	};

	if (_BEMDomain[domain]) {
		if (matrices & Matrices::M) {
			ESINFO(ERROR) << "BEM not support computation of matrix M.";
		}
		if (matrices & Matrices::R) {
			ESINFO(ERROR) << "BEM not support computation of matrix R.";
		}
		processBEM(domain, matrices);
		_instance->K[domain].ConvertDenseToCSR(1);
		_instance->K[domain].RemoveLower();
		boundary(matrices & Matrices::f);
	} else {
		auto enodes = _mesh->elements->procNodes->cbegin() + _mesh->elements->elementsDistribution[domain];
		for (eslocal e = _mesh->elements->elementsDistribution[domain]; e < _mesh->elements->elementsDistribution[domain + 1]; ++e, ++enodes) {
			processElement(domain, matrices, e, Ke, Me, Re, fe);
			switch (getMatrixType(domain)) {
			case MatrixType::REAL_UNSYMMETRIC: fullInsert(enodes, 1); break;
			default: upperInsert(enodes, 1);
			}
		}

		boundary(matrices);

		_instance->K[domain].mtype = getMatrixType(domain);
		_instance->M[domain].mtype = MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;
	}
}

void Physics::makeStiffnessMatricesRegular(FETI_REGULARIZATION regularization, size_t scSize, bool ortogonalCluster)
{
	#pragma omp parallel for
	for (size_t d = 0; d < _instance->domains; d++) {
		makeStiffnessMatrixRegular(regularization, scSize, d, ortogonalCluster);
		ESINFO(PROGRESS3) << Info::plain() << ".";
	}
	ESINFO(PROGRESS3);
}

void Physics::makeStiffnessMatrixRegular(FETI_REGULARIZATION regularization, size_t scSize, size_t domain, bool ortogonalCluster)
{
	switch (regularization) {

	case FETI_REGULARIZATION::ANALYTIC:
		analyticRegularization(domain, ortogonalCluster);
		_instance->RegMat[domain].RemoveLower();
		_instance->K[domain].MatAddInPlace(_instance->RegMat[domain], 'N', 1);
		_instance->RegMat[domain].ConvertToCOO(1);
		break;

	case FETI_REGULARIZATION::ALGEBRAIC:
		switch (_instance->K[domain].mtype) {
			double norm;
			eslocal defect;

		case MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE:
			_instance->K[domain].get_kernel_from_K(_instance->K[domain], _instance->RegMat[domain], _instance->N1[domain], norm, defect, domain, scSize);
			break;

		case MatrixType::REAL_UNSYMMETRIC:
			_instance->K[domain].get_kernels_from_nonsym_K(_instance->K[domain], _instance->RegMat[domain], _instance->N1[domain], _instance->N2[domain], norm, defect, domain, scSize);
			break;

		default:
			ESINFO(ERROR) << "Unknown matrix type for regularization.";
		}
		break;
	}
}

double Physics::sumSquares(const std::vector<std::vector<double> > &data, SumRestriction restriction) const
{
	ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: implement sumSquares.";
	return 0;
}

void Physics::assembleB1(bool withRedundantMultipliers, bool withGluing, bool withScaling)
{
	_constraints->B1DirichletInsert(*_step);
	if (withGluing) {
		_constraints->B1GlueElements();
	}
}

void Physics::updateDirichletInB1(bool withRedundantMultipliers)
{
	_constraints->B1DirichletUpdate(*_step);
}

void Physics::updateDuplicity()
{
	_constraints->B1DuplicityUpdate();
}

void Physics::assembleB0FromCorners()
{
	_constraints->B0Corners();
}

void Physics::assembleB0FromKernels(const std::vector<SparseMatrix> &kernels)
{
	_constraints->B0Kernels(kernels);
}

void Physics::smoothstep(double &smoothStep, double &derivation, double edge0, double edge1, double value, size_t order)
{
	value = std::max(0.0, std::min((value - edge0) / (edge1 - edge0), 1.0));

	switch (order) {
	case 0:
		smoothStep = value;
		if (value == 0 || value == 1) {
			derivation = 0;
		} else {
			derivation = 1;
		}
		break;
	case 1:
		smoothStep = -2 * pow(value, 3) + 3 * pow(value, 2);
		derivation = -6 * pow(value, 2) + 6 * value;
		break;
	case 2:
		smoothStep = 6 * pow(value, 5) - 15 * pow(value, 4) + 10 * pow(value, 3);
		derivation = 30 * pow(value, 4) - 60*  pow(value, 3) + 30 * pow(value, 2);
		break;
	case 3:
		smoothStep = -20 * pow(value, 7) + 70 * pow(value, 6) - 84 * pow(value, 5) + 35 * pow(value, 4);
		derivation = -140 * pow(value, 6) + 420 * pow(value, 5) - 420 * pow(value, 4) + 140 * pow(value, 3);
		break;
	case 4:
		smoothStep = 70 * pow(value, 9) - 315 * pow(value, 8) + 540 * pow(value, 7) - 420 * pow(value, 6) + 126 * pow(value, 5);
		derivation = 630 * pow(value, 8) - 2520 * pow(value, 7) + 3780 * pow(value, 6) - 2520 * pow(value, 5) + 630 * pow(value, 4);
		break;
	case 5:
		smoothStep = -252 * pow(value, 11) + 1386 * pow(value, 10) - 3080 * pow(value, 9) + 3465 * pow(value, 8) - 1980 * pow(value, 7) + 462 * pow(value, 6);
		derivation = -2772 * pow(value, 10) + 13860 * pow(value, 9) - 27720 * pow(value, 8) + 27720 * pow(value, 7) - 13860 * pow(value, 6) + 2772 * pow(value, 5);
		break;
	case 6:
		smoothStep = 924 * pow(value, 13) - 6006 * pow(value, 12) + 16380 * pow(value, 11) - 24024 * pow(value, 10) + 20020 * pow(value, 9) - 9009 * pow(value, 8) + 1716 * pow(value, 7);
		derivation = 12012 * pow(value, 12) - 72072 * pow(value, 11) + 180180 * pow(value, 10) - 240240 * pow(value, 9) + 180180 * pow(value, 8) - 72072 * pow(value, 7) + 12012 * pow(value, 6);
		break;
	}

	derivation /= edge1 - edge0;
}


