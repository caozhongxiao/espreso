
#include "loader.h"
#include "workbench/workbench.h"

#include "../basis/containers/point.h"
#include "../basis/containers/serializededata.h"
#include "../basis/utilities/utils.h"
#include "../basis/utilities/communication.h"
#include "../basis/logging/timeeval.h"

#include "../config/ecf/root.h"

#include "../mesh/mesh.h"
#include "../mesh/preprocessing/meshpreprocessing.h"
#include "../mesh/elements/element.h"
#include "../mesh/store/nodestore.h"
#include "../mesh/store/elementstore.h"
#include "../mesh/store/elementsregionstore.h"
#include "../mesh/store/boundaryregionstore.h"
#include "../old/input/loader.h"

#include "../wrappers/metis/wparmetis.h"

#include <numeric>
#include <algorithm>

#include <fstream>

using namespace espreso;

void BalancedLoader::load(const ECFRoot &configuration, Mesh &mesh, int MPIrank, int MPIsize)
{
	switch (configuration.input) {
	case INPUT_FORMAT::WORKBENCH:
		WorkbenchLoader::load(configuration, mesh);
		mesh.update();
		break;
	default:
		input::OldLoader::load(configuration, *mesh.mesh, MPIrank, MPIsize);
		mesh.load();
		break;
	}
}

void BalancedLoader::loadDistributedMesh(const ECFRoot &configuration, DistributedMesh &dMesh, Mesh &mesh)
{
	BalancedLoader(configuration, dMesh, mesh);
}

BalancedLoader::BalancedLoader(const ECFRoot &configuration, DistributedMesh &dMesh, Mesh &mesh)
: _configuration(configuration), _dMesh(dMesh), _mesh(mesh)
{
	ESINFO(OVERVIEW) << "Balance distributed mesh.";
	TimeEval timing("Load distributed mesh");
	timing.totalTime.startWithBarrier();

	TimeEvent tdistribution("distribute mesh across processes"); tdistribution.start();
	distributeMesh();
	tdistribution.end(); timing.addEvent(tdistribution);
	ESINFO(PROGRESS2) << "Balanced loader:: data balanced.";

	TimeEvent tesort("sort elements accordint to the first node"); tesort.start();
//	sortElementsVariousTargets();
	sortElementsManual();
	SFC();
	tesort.end(); timing.addEvent(tesort);
	ESINFO(PROGRESS2) << "Balanced loader:: elements sorted.";

	TimeEvent telements("fill elements"); telements.start();
	fillElements();
	telements.end(); timing.addEvent(telements);
	ESINFO(PROGRESS2) << "Balanced loader:: elements filled.";

	TimeEvent tcoordinates("fill coordinates"); tcoordinates.start();
	fillCoordinates();
	tcoordinates.end(); timing.addEvent(tcoordinates);
	ESINFO(PROGRESS2) << "Balanced loader:: coordinates filled.";

//	TimeEvent tpolish("polishing of decomposition"); tpolish.start();
//	polishSFC();
//	tpolish.end(); timing.addEvent(tpolish);
//	ESINFO(PROGRESS2) << "Balanced loader:: decomposition polished.";

	timing.totalTime.endWithBarrier();
	timing.printStatsMPI();
}

void BalancedLoader::distributeMesh()
{
	// DISTRIBUTE NODES
	eslocal myMaxID = 0, maxID;
	int sorted = std::is_sorted(_dMesh.nIDs.begin(), _dMesh.nIDs.end()), allSorted;
	std::vector<eslocal> permutation;
	if (_dMesh.nIDs.size()) {
		if (sorted) {
			myMaxID = _dMesh.nIDs.back();
		} else {
			permutation.resize(_dMesh.nIDs.size());
			std::iota(permutation.begin(), permutation.end(), 0);
			std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) { return _dMesh.nIDs[i] < _dMesh.nIDs[j]; });
			myMaxID = _dMesh.nIDs[permutation.back()];
		}
	}

	MPI_Allreduce(&myMaxID, &maxID, sizeof(eslocal), MPI_BYTE, MPITools::operations().max, environment->MPICommunicator);
	MPI_Allreduce(&sorted, &allSorted, 1, MPI_INT, MPI_MIN, environment->MPICommunicator);

	if (allSorted) {
		if (environment->MPIsize == 1) {
			_nDistribution = { 0, _dMesh.nIDs.size() };
		} else {
			std::vector<size_t> cCurrent = Communication::getDistribution(_dMesh.nIDs.size(), MPITools::operations().sizeToOffsetsSize_t);
			_nDistribution = tarray<eslocal>::distribute(environment->MPIsize, cCurrent.back());

			if (!Communication::balance(_dMesh.nIDs, cCurrent, _nDistribution)) {
				ESINFO(ERROR) << "ESPRESO internal error: balance node IDs.";
			}
			if (!Communication::balance(_dMesh.coordinates, cCurrent, _nDistribution)) {
				ESINFO(ERROR) << "ESPRESO internal error: balance coordinates.";
			}
		}
	} else {
		if (environment->MPIsize == 1) {
			_nDistribution = { 0, _dMesh.nIDs.size() };
			std::vector<eslocal> nIDs;
			std::vector<Point> nPoints;
			nIDs.reserve(_dMesh.nIDs.size());
			nPoints.reserve(_dMesh.nIDs.size());
			for (size_t i = 0; i < permutation.size(); i++) {
				nIDs.push_back(_dMesh.nIDs[permutation[i]]);
				nPoints.push_back(_dMesh.coordinates[permutation[i]]);
			}
			_dMesh.nIDs.swap(nIDs);
			_dMesh.coordinates.swap(nPoints);
		} else {
			_nDistribution = tarray<eslocal>::distribute(environment->MPIsize, maxID + 1);
			std::vector<std::vector<eslocal> > sIDs, rIDs;
			std::vector<std::vector<Point> > sCoordinates, rCoordinates;
			std::vector<int> targets;
			for (int r = 0; r < environment->MPIsize; r++) {
				auto begin = std::lower_bound(permutation.begin(), permutation.end(), _nDistribution[r], [&] (eslocal i, const size_t &ID) { return _dMesh.nIDs[i] < ID; });
				auto end = std::lower_bound(permutation.begin(), permutation.end(), _nDistribution[r + 1], [&] (eslocal i, const size_t &ID) { return _dMesh.nIDs[i] < ID; });
				if (begin != end) {
					sIDs.push_back({});
					sCoordinates.push_back({});
					targets.push_back(r);
				}
				for (size_t n = begin - permutation.begin(); n < end - permutation.begin(); ++n) {
					sIDs.back().push_back(_dMesh.nIDs[permutation[n]]);
					sCoordinates.back().push_back(_dMesh.coordinates[permutation[n]]);
				}
			}

			if (!Communication::sendVariousTargets(sIDs, rIDs, targets)) {
				ESINFO(ERROR) << "ESPRESO internal error: distribute not sorted node IDs.";
			}
			if (!Communication::sendVariousTargets(sCoordinates, rCoordinates, targets)) {
				ESINFO(ERROR) << "ESPRESO internal error: distribute not sorted node coordinates.";
			}

			for (size_t r = 1; r < rIDs.size(); r++) {
				rIDs[0].insert(rIDs[0].end(), rIDs[r].begin(), rIDs[r].end());
				rCoordinates[0].insert(rCoordinates[0].end(), rCoordinates[r].begin(), rCoordinates[r].end());
			}

			permutation.resize(rIDs[0].size());
			std::iota(permutation.begin(), permutation.end(), 0);
			std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) { return rIDs[0][i] < rIDs[0][j]; });
			_dMesh.nIDs.clear();
			_dMesh.coordinates.clear();
			_dMesh.nIDs.reserve(permutation.size());
			_dMesh.coordinates.reserve(permutation.size());
			for (size_t n = 0; n < permutation.size(); n++) {
				_dMesh.nIDs.push_back(rIDs[0][permutation[n]]);
				_dMesh.coordinates.push_back(rCoordinates[0][permutation[n]]);
			}
		}
	}

	// DISTRIBUTE ELEMENTS

	myMaxID = 0;
	sorted = std::is_sorted(_dMesh.edata.begin(), _dMesh.edata.end(), [] (const EData &e1, const EData &e2) { return e1.id < e2.id; }), allSorted;
	if (_dMesh.edata.size()) {
		if (sorted) {
			myMaxID = _dMesh.edata.back().id;
		} else {
			permutation.resize(_dMesh.edata.size());
			std::iota(permutation.begin(), permutation.end(), 0);
			std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) { return _dMesh.edata[i].id < _dMesh.edata[j].id; });
			myMaxID = _dMesh.edata[permutation.back()].id;
		}
	}

	MPI_Allreduce(&myMaxID, &maxID, sizeof(eslocal), MPI_BYTE, MPITools::operations().max, environment->MPICommunicator);
	MPI_Allreduce(&sorted, &allSorted, 1, MPI_INT, MPI_MIN, environment->MPICommunicator);

	if (allSorted) {
		std::vector<size_t> eCurrent = Communication::getDistribution(_dMesh.esize.size(), MPITools::operations().sizeToOffsetsSize_t);
		std::vector<size_t> eTarget = tarray<eslocal>::distribute(environment->MPIsize, eCurrent.back());

		std::vector<size_t> nCurrent = Communication::getDistribution(_dMesh.enodes.size(), MPITools::operations().sizeToOffsetsSize_t);
		std::vector<size_t> nTarget;

		if (environment->MPIrank == 0) {
			nTarget.push_back(0);
		}

		size_t nodeOffset = nCurrent[environment->MPIrank];
		size_t eTargetIndex = std::lower_bound(eTarget.begin(), eTarget.end(), eCurrent[environment->MPIrank] + 1) - eTarget.begin();
		for (size_t n = 0; n < _dMesh.esize.size(); ++n) {
			nodeOffset += _dMesh.esize[n];
			if (eCurrent[environment->MPIrank] + n + 1 == eTarget[eTargetIndex]) {
				nTarget.push_back(nodeOffset);
				++eTargetIndex;
			}
		}
		Communication::allGatherUnknownSize(nTarget);

		if (!Communication::balance(_dMesh.esize, eCurrent, eTarget)) {
			ESINFO(ERROR) << "ESPRESO internal error: balance element sizes.";
		}
		if (!Communication::balance(_dMesh.edata, eCurrent, eTarget)) {
			ESINFO(ERROR) << "ESPRESO internal error: balance element data.";
		}
		if (!Communication::balance(_dMesh.enodes, nCurrent, nTarget)) {
			ESINFO(ERROR) << "ESPRESO internal error: balance element nodes.";
		}

		_eDistribution = eTarget;
	} else {
		if (environment->MPIsize == 1) {
			_eDistribution = { 0, _dMesh.esize.size() };
			std::vector<eslocal> eSize, eNodes;
			std::vector<EData> eData;
			eSize.reserve(_dMesh.esize.size());
			eData.reserve(_dMesh.esize.size());
			eNodes.reserve(_dMesh.enodes.size());
			std::vector<eslocal> edist = { 0 };
			edist.reserve(_dMesh.esize.size() + 1);
			for (size_t e = 0; e < _dMesh.esize.size(); e++) {
				edist.push_back(edist.back() + _dMesh.esize[e]);
			}

			for (size_t i = 0; i < permutation.size(); i++) {
				eData.push_back(_dMesh.edata[permutation[i]]);
				eSize.push_back(_dMesh.esize[permutation[i]]);
				eNodes.insert(eNodes.end(), _dMesh.enodes.begin() + edist[permutation[i]], _dMesh.enodes.begin() + edist[permutation[i] + 1]);
			}

			_dMesh.esize.swap(eSize);
			_dMesh.edata.swap(eData);
			_dMesh.enodes.swap(eNodes);
		} else {
			_eDistribution = tarray<eslocal>::distribute(environment->MPIsize, maxID + 1);
			std::vector<std::vector<eslocal> > sSize, sNodes, rSize, rNodes;
			std::vector<std::vector<EData> > sEData, rEData;
			std::vector<int> targets;
			std::vector<eslocal> edist = { 0 };
			edist.reserve(_dMesh.esize.size() + 1);
			for (size_t e = 0; e < _dMesh.esize.size(); e++) {
				edist.push_back(edist.back() + _dMesh.esize[e]);
			}

			if (sorted) {
				for (int r = 0; r < environment->MPIsize; r++) {
					auto begin = std::lower_bound(_dMesh.edata.begin(), _dMesh.edata.end(), _eDistribution[r], [&] (EData &edata, const size_t &ID) { return edata.id < ID; });
					auto end = std::lower_bound(_dMesh.edata.begin(), _dMesh.edata.end(), _eDistribution[r + 1], [&] (EData &edata, const size_t &ID) { return edata.id < ID; });
					if (begin != end) {
						sSize.push_back({});
						sNodes.push_back({});
						sEData.push_back({});
						targets.push_back(r);
						size_t b = begin - _dMesh.edata.begin();
						size_t e = end - _dMesh.edata.begin();
						sSize.back().insert(sSize.back().end(), _dMesh.esize.begin() + b, _dMesh.esize.begin() + e);
						sEData.back().insert(sEData.back().end(), _dMesh.edata.begin() + b, _dMesh.edata.begin() + e);
						sNodes.back().insert(sNodes.back().end(), _dMesh.enodes.begin() + edist[b], _dMesh.enodes.begin() + edist[e]);
					}
				}
			} else {
				for (int r = 0; r < environment->MPIsize; r++) {
					auto begin = std::lower_bound(permutation.begin(), permutation.end(), _eDistribution[r], [&] (eslocal i, const size_t &ID) { return _dMesh.edata[i].id < ID; });
					auto end = std::lower_bound(permutation.begin(), permutation.end(), _eDistribution[r + 1], [&] (eslocal i, const size_t &ID) { return _dMesh.edata[i].id < ID; });
					if (begin != end) {
						sSize.push_back({});
						sNodes.push_back({});
						sEData.push_back({});
						targets.push_back(r);
					}
					for (size_t n = begin - permutation.begin(); n < end - permutation.begin(); ++n) {
						sSize.back().push_back(_dMesh.esize[permutation[n]]);
						sEData.back().push_back(_dMesh.edata[permutation[n]]);
						sNodes.back().insert(sNodes.back().end(), _dMesh.enodes.begin() + edist[permutation[n]], _dMesh.enodes.begin() + edist[permutation[n] + 1]);
					}
				}
			}

			if (!Communication::sendVariousTargets(sSize, rSize, targets)) {
				ESINFO(ERROR) << "ESPRESO internal error: distribute not sorted elements sizes.";
			}
			if (!Communication::sendVariousTargets(sEData, rEData, targets)) {
				ESINFO(ERROR) << "ESPRESO internal error: distribute not sorted element data.";
			}
			if (!Communication::sendVariousTargets(sNodes, rNodes, targets)) {
				ESINFO(ERROR) << "ESPRESO internal error: distribute not sorted element nodes.";
			}

			for (size_t r = 1; r < rSize.size(); r++) {
				rSize[0].insert(rSize[0].end(), rSize[r].begin(), rSize[r].end());
				rEData[0].insert(rEData[0].end(), rEData[r].begin(), rEData[r].end());
				rNodes[0].insert(rNodes[0].end(), rNodes[r].begin(), rNodes[r].end());
			}

			permutation.resize(rSize[0].size());
			std::iota(permutation.begin(), permutation.end(), 0);
			std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) { return rEData[0][i].id < rEData[0][j].id; });

			edist = std::vector<eslocal>({ 0 });
			edist.reserve(rSize[0].size() + 1);
			for (size_t e = 0; e < rSize[0].size(); e++) {
				edist.push_back(edist.back() + rSize[0][e]);
			}

			_dMesh.esize.clear();
			_dMesh.enodes.clear();
			_dMesh.edata.clear();
			_dMesh.esize.reserve(permutation.size());
			_dMesh.edata.reserve(permutation.size());
			_dMesh.enodes.reserve(rNodes[0].size());
			for (size_t n = 0; n < permutation.size(); n++) {
				_dMesh.esize.push_back(rSize[0][permutation[n]]);
				_dMesh.edata.push_back(rEData[0][permutation[n]]);
				_dMesh.enodes.insert(_dMesh.enodes.end(), rNodes[0].begin() + edist[permutation[n]], rNodes[0].begin() + edist[permutation[n] + 1]);
			}
		}
	}
}

static void rotateD2(size_t n, size_t &x, size_t &y, int rx, int ry) {
	if (ry == 0) {
		if (rx == 1) {
			x = n - 1 - x;
			y = n - 1 - y;
		}
		std::swap(x, y);
	}
}

static void rotateD3(size_t n, size_t &x, size_t &y, size_t &z, int rx, int ry, int rz) {
	if (rz == 0) {
		if (rx == 1) {
			std::swap(x, z);
			x = n - 1 - x;
			z = n - 1 - z;
		} else {
			if (ry == 1) {
				y = n - 1 - y;
				z = n - 1 - z;
			}
			std::swap(y, z);
		}
	} else {
		std::swap(x, y);
		if (ry == 1) {
			x = n - 1 - x;
			y = n - 1 - y;
		}
	}
}

static double D2toD1(size_t n, size_t x, size_t y) {
	int rx, ry;
	double d = 0;
	for (size_t s = n / 2, depth = 2; s > 0; s /= 2, depth = depth << 1) {
		rx = (x & s) > 0;
		ry = (y & s) > 0;
		d += 1. / (depth * depth) * ((3 * rx) ^ ry);

		rotateD2(n, x, y, rx, ry);
	}
	return d;
};

static double D3toD1(size_t n, size_t x, size_t y, size_t z) {
	int rx, ry, rz;
	double d = 0;

	for (size_t s = n / 2, depth = 2; s > 0; s /= 2, depth = depth << 1) {
		rx = (x & s) > 0;
		ry = (y & s) > 0;
		rz = (z & s) > 0;
		d += 1. / (depth * depth * depth) * (4 * ry + (3 * (ry ^ rx)) ^ rz);

		rotateD3(n, x, y, z, rx, ry, rz);
	}
	return d;
};

static void D1toD2(size_t n, size_t d, size_t &x, size_t &y) {
	int rx, ry;
	x = y = 0;
	for (size_t s = 1, t = d; s < n; s *= 2, t /= 4) {
		rx = 1 & (t / 2);
		ry = 1 & (t ^ rx);

		rotateD2(s, x, y, rx, ry);

		x += s * rx;
		y += s * ry;
	}
};

static void D1toD3(size_t n, size_t d, size_t &x, size_t &y, size_t &z) {
	int rx, ry, rz;
	x = y = z = 0;
	for (size_t s = 1, t = d; s < n; s *= 2, t /= 8) {
		rx = (1 & (t / 4)) ^ (1 & (t / 2));
		ry = (1 & (t / 4));
		rz = (1 & (t / 2)) ^ (1 & t);

		rotateD3(s, x, y, z, rx, ry, rz);

		x += s * rx;
		y += s * ry;
		z += s * rz;
	}
};


void BalancedLoader::SFC()
{
	Point min = _dMesh.coordinates.front(), max = _dMesh.coordinates.front();

	for (size_t i = 1; i < _dMesh.coordinates.size(); i++) {
		min.x = std::min(min.x, _dMesh.coordinates[i].x);
		min.y = std::min(min.y, _dMesh.coordinates[i].y);
		min.z = std::min(min.z, _dMesh.coordinates[i].z);
		max.x = std::max(max.x, _dMesh.coordinates[i].x);
		max.y = std::max(max.y, _dMesh.coordinates[i].y);
		max.z = std::max(max.z, _dMesh.coordinates[i].z);
	}

	Point origin, size;
	MPI_Allreduce(&min.x, &origin.x, 3, MPI_DOUBLE, MPI_MIN, environment->MPICommunicator);
	MPI_Allreduce(&max.x, &size.x, 3, MPI_DOUBLE, MPI_MAX, environment->MPICommunicator);

	_mesh.nodes->min = origin;
	_mesh.nodes->max = size;

	size -= origin - Point(1e-6, 1e-6, 1e-6);

	eslocal dimension = 0;
	switch (_configuration.physics) {
	case PHYSICS::HEAT_TRANSFER_2D:
	case PHYSICS::STRUCTURAL_MECHANICS_2D:
	case PHYSICS::SHALLOW_WATER_2D:
		dimension = 2;
		break;
	case PHYSICS::HEAT_TRANSFER_3D:
	case PHYSICS::STRUCTURAL_MECHANICS_3D:
		dimension = 3;
		break;
	}

	size_t N = 1, NN;
	while (pow(N, dimension) < environment->MPIsize) {
		N = N << 1;
	}
	N = 2;
	NN = N;

	std::vector<size_t> sumoffset, nextsumoffset = { 0 };
	std::vector<std::vector<size_t> > refined = { {0} }, refinedxyz = { {0} };
	std::vector<eslocal> scounts(pow(N, dimension)), rcounts(pow(N, dimension));
	std::vector<double> partition(_dMesh.esize.size());
	std::vector<size_t> asum(pow(N, dimension) + 1), ideal = tarray<size_t>::distribute(environment->MPIsize, _eDistribution.back());
	std::vector<eslocal> permutation(partition.size());

	_sfcbounds.resize(environment->MPIsize + 1);
	_sfcbounds.back() = 1;
	std::iota(permutation.begin(), permutation.end(), 0);

	std::vector<eslocal> edist({ 0 });
	edist.reserve(_dMesh.esize.size() + 1);
	for (size_t e = 0; e < _dMesh.esize.size(); e++) {
		edist.push_back(edist.back() + _dMesh.esize[e]);
	}

	double PRECISION = 0.02 * std::log2(environment->MPIsize);

	std::vector<std::pair<size_t, size_t> > sfcboundary(environment->MPIsize + 1);
	sfcboundary.back().first = sfcboundary.front().first = N;
	sfcboundary.back().second = pow(N, dimension) - 1;
	while (refined.back().size()) {
		sumoffset.swap(nextsumoffset);
		nextsumoffset.clear();
		refined.push_back({});

		auto begin = permutation.begin();
		auto end = permutation.begin();
		for (auto q = refined[refined.size() - 2].begin(), offset = sumoffset.begin(); q != refined[refined.size() - 2].end(); ++q, ++offset) {
			std::fill(scounts.begin(), scounts.end(), 0);

			double lower = *q * 1. / pow(N / 2, dimension);
			double upper = lower + 1. / pow(N / 2, dimension);
			if (N == NN) {
				lower = 0;
				upper = 1;
			}

			begin = std::lower_bound(end, permutation.end(), lower, [&] (eslocal i, double b) { return partition[i] < b; });
			end = std::lower_bound(begin, permutation.end(), upper, [&] (eslocal i, double b) { return partition[i] < b; });

			if (dimension == 2) {
				for (auto i = begin; i != end; ++i) {
					size_t x = std::floor(N * (_dMesh.coordinates[_dMesh.enodes[edist[*i]] - _nDistribution[environment->MPIrank]].x - origin.x) / size.x);
					size_t y = std::floor(N * (_dMesh.coordinates[_dMesh.enodes[edist[*i]] - _nDistribution[environment->MPIrank]].y - origin.y) / size.y);
					++scounts[((partition[*i] = D2toD1(N, x, y)) - lower) / (1. / pow(N, dimension))];
				}
			}

			if (dimension == 3) {
				for (auto i = begin; i != end; ++i) {
					const Point &p = _dMesh.coordinates[_dMesh.enodes[edist[*i]] - _nDistribution[environment->MPIrank]];
					size_t x = std::floor(N * (p.x - origin.x) / size.x);
					size_t y = std::floor(N * (p.y - origin.y) / size.y);
					size_t z = std::floor(N * (p.z - origin.z) / size.z);
					++scounts[((partition[*i] = D3toD1(N, x, y, z)) - lower) / (1. / pow(N, dimension))];
				}
			}

			std::sort(begin, end, [&] (eslocal i, eslocal j) { return partition[i] < partition[j]; });

			MPI_Allreduce(scounts.data(), rcounts.data(), sizeof(eslocal) * scounts.size(), MPI_BYTE, MPITools::operations().sum, environment->MPICommunicator);

			asum[0] = *offset;
			for (size_t i = 0; i < rcounts.size(); i++) {
				asum[i + 1] = asum[i] + rcounts[i];
			}

			auto ibegin = std::lower_bound(ideal.begin(), ideal.end(), asum.front());
			auto iend = std::lower_bound(ideal.begin(), ideal.end(), asum.back());
			for (auto i = ibegin; i != iend; ++i) {
				auto up = std::lower_bound(asum.begin(), asum.end(), *i);
					if (up != asum.begin()) {
					auto bottom = up - 1;
					double scale = (double)(*i - *bottom) / (*up - *bottom);
					_sfcbounds[i - ideal.begin()] = 1. / pow(N, dimension) * ((*q * pow(2, dimension)) + bottom - asum.begin()) + scale * 1. / pow(N, dimension);
					sfcboundary[i - ideal.begin()].first = N;
					sfcboundary[i - ideal.begin()].second = (*q * pow(2, dimension)) + bottom - asum.begin();
					if (*up - *bottom > PRECISION * (_eDistribution.back() / environment->MPIsize)) {
						if (refined.back().size() == 0 || refined.back().back() != (*q * pow(2, dimension)) + bottom - asum.begin()) {
							refined.back().push_back((*q * pow(2, dimension)) + bottom - asum.begin());
							nextsumoffset.push_back(*bottom);
						}
					}
					if ((*q * pow(2, dimension)) + bottom - asum.begin() == pow(N, dimension) - 1) {
						sfcboundary.back().first = N << 1;
						sfcboundary.back().second = pow(N << 1, dimension) - 1;
					}
				}
			}
		}
		if (refined.back().size()) {
			N = N << 1;
		}
		scounts.resize(pow(2, dimension));
		rcounts.resize(pow(2, dimension));
		asum.resize(pow(2, dimension) + 1);
	}

	size_t x, y, z = 0, xx, yy, zz, depth;
	for (size_t i = 1; i < refined.size(); i++) {
		refinedxyz.push_back({});
		for (size_t j = 0; j < refined[i].size(); j++) {
			if (dimension == 2) {
				D1toD2(NN << (i - 1), refined[i][j], x, y);
				refinedxyz[i].push_back(y * (NN << (i - 1)) + x);
			}
			if (dimension == 3) {
				D1toD3(NN << (i - 1), refined[i][j], x, y, z);
				refinedxyz[i].push_back(z * (NN << (i - 1)) * (NN << (i - 1)) + y * (NN << (i - 1)) + x);
			}
		}
		std::sort(refinedxyz[i].begin(), refinedxyz[i].end());
	}

	if (environment->MPIrank == 0) {
		std::cout << "NN: " << NN << ", N: " << N << "\n";
		for (size_t i = 0; i < refined.size(); i++) {
			std::cout << refined[i];
		}
//		for (size_t i = 0; i < refinedxyz.size(); i++) {
//			std::cout << refinedxyz[i];
//		}
//		std::cout << _sfcbounds;
//
//		for (size_t i = 0; i < sfcboundary.size(); i++) {
//			D1toD2(sfcboundary[i].first, sfcboundary[i].second, x, y);
////			std::cout << sfcboundary[i].first << "::" << sfcboundary[i].second << " -> [" << x << "," << y << "]\n";
//		}
	}

	if (environment->MPIrank == 0) {
		std::ofstream os("SFC.vtk");
		os << "# vtk DataFile Version 2.0\n";
		os << "EXAMPLE\n";
		os << "ASCII\n";
		os << "DATASET UNSTRUCTURED_GRID\n\n";

		os << "POINTS " << pow(N + 1, dimension) << " float\n";

		for (size_t k = 0; k <= N; k++) {
			for (size_t j = 0; j <= N; j++) {
				for (size_t i = 0; i <= N; i++) {
					os << i * size.x / N << " " << j * size.y / N << " " << k * size.z / N << " \n";
				}
			}
		}
		os << "\n";

		std::vector<size_t> ccount = { 0 };
		size_t cells = 0, prev = pow(NN, dimension);
		for (size_t d = 1; d < refinedxyz.size(); d++) {
			std::cout << "prev: " << prev << ", xy: " << refinedxyz[d].size() << "\n";
			cells += prev - refinedxyz[d].size();
			prev = pow(2, dimension) * refinedxyz[d].size();
		}
		os << "CELLS " << cells << " " << cells + pow(2, dimension) * cells << "\n";

		std::cout << sfcboundary;
		size_t n = sfcboundary[environment->MPIrank].first;
		size_t cell = sfcboundary[environment->MPIrank].second;
		size_t level = 1;
		while ((NN << level) <= n) ++level;
		std::vector<std::vector<size_t>::const_iterator> its = { refined.front().begin() };

//		std::cout << "n: " << n << ", cell: " << cell << ", level: " << level << "\n";
		cell *= pow(N / n, dimension);
		for (size_t i = 1; i < refined.size(); i++) {
			its.push_back(std::lower_bound(refined[i].begin(), refined[i].end(), cell / (size_t)pow(1 << (refined.size() - i - 1), dimension)));
		}
		cell = sfcboundary[environment->MPIrank].second;

		size_t index = 0;
		for (int r = 0; r < environment->MPIsize; r++) {
			while (n != sfcboundary[r + 1].first || cell <= sfcboundary[r + 1].second) {
				while (its[level] != refined[level].end() && *its[level] == cell) {
					++level;
					n = n << 1;
					cell *= pow(2, dimension);
				}

				while (n > NN && its[level - 1] != refined[level - 1].end() && *its[level - 1] < cell / (size_t)pow(2, dimension)) {
					++its[level - 1];
					if (its[level - 1] == refined[level - 1].end() || *its[level - 1] != cell / (size_t)pow(2, dimension)) {
						--level;
						n = n >> 1;
						cell /= pow(2, dimension);
					}
				}

				if (dimension == 2) {
					D1toD2(n, cell, x, y);
				}
				if (dimension == 3) {
					D1toD3(n, cell, x, y, z);
				}

//				std::cout << "level: " << level << ", n: " << n << ", cell: " << cell << " - " << z * n * n + y * n + x << "\n";

				size_t row = N + 1;
//				std::cout << "ROW: " << row << "\n";
				if (dimension == 2) {
					os << "4 ";
					os << row * (N / n) * (y + 0) + (N / n) * (x + 0) << " ";
					os << row * (N / n) * (y + 0) + (N / n) * (x + 1) << " ";
					os << row * (N / n) * (y + 1) + (N / n) * (x + 1) << " ";
					os << row * (N / n) * (y + 1) + (N / n) * (x + 0) << "\n";
				}

				if (dimension == 3) {
//					std::cout << "x: " << x << ", y: " << y << ", z: " << z << ", N / n: " << N / n << "\n";
					os << "8 ";
					os << row * row * (N / n) * (z + 0) + row * (N / n) * (y + 0) + (N / n) * (x + 0) << " ";
					os << row * row * (N / n) * (z + 0) + row * (N / n) * (y + 0) + (N / n) * (x + 1) << " ";
					os << row * row * (N / n) * (z + 0) + row * (N / n) * (y + 1) + (N / n) * (x + 1) << " ";
					os << row * row * (N / n) * (z + 0) + row * (N / n) * (y + 1) + (N / n) * (x + 0) << " ";
					os << row * row * (N / n) * (z + 1) + row * (N / n) * (y + 0) + (N / n) * (x + 0) << " ";
					os << row * row * (N / n) * (z + 1) + row * (N / n) * (y + 0) + (N / n) * (x + 1) << " ";
					os << row * row * (N / n) * (z + 1) + row * (N / n) * (y + 1) + (N / n) * (x + 1) << " ";
					os << row * row * (N / n) * (z + 1) + row * (N / n) * (y + 1) + (N / n) * (x + 0) << "\n";
				}

				++cell;
				++index;
			}
			ccount.push_back(index);
		}
		std::cout << ccount;

		os << "\n";

		os << "CELL_TYPES " << cells << "\n";
		for (size_t i = 0; i < cells; i++) {
			if (dimension == 2) {
				os << "9\n";
			}
			if (dimension == 3) {
				os << "12\n";
			}
		}
		os << "\n";

		os << "CELL_DATA " << cells << "\n";
		os << "SCALARS MPI int 1\n";
		os << "LOOKUP_TABLE default\n";
		for (int r = 0; r < environment->MPIsize; r++) {
			for (size_t i = ccount[r]; i < ccount[r + 1]; i++) {
				os << r << "\n";
			}
		}
		os << "\n";

		os.close();
	}

	std::vector<std::pair<size_t, size_t> > neighbors, potential;
	std::vector<std::pair<double, double> > intervals;

	Communication::serialize([&] () {
//		std::cout << " --- " << environment->MPIrank << " -- \n";
		size_t n = sfcboundary[environment->MPIrank].first;
		size_t cell = sfcboundary[environment->MPIrank].second;
		size_t level = 1;
		while ((NN << level) <= n) ++level;
		std::vector<std::vector<size_t>::const_iterator> its = { refined.front().begin() };

//		std::cout << "n: " << n << ", cell: " << cell << ", level: " << level << "\n";
		cell *= pow(N / n, dimension);
		for (size_t i = 1; i < refined.size(); i++) {
			its.push_back(std::lower_bound(refined[i].begin(), refined[i].end(), cell / (size_t)pow(1 << (refined.size() - i - 1), dimension)));
		}
		cell = sfcboundary[environment->MPIrank].second;
		while (n != sfcboundary[environment->MPIrank + 1].first || cell <= sfcboundary[environment->MPIrank + 1].second) {
			while (its[level] != refined[level].end() && *its[level] == cell) {
				++level;
				n = n << 1;
				cell *= pow(2, dimension);
			}

			while (n > NN && its[level - 1] != refined[level - 1].end() && *its[level - 1] < cell / (size_t)pow(2, dimension)) {
//				std::cout << "DOWN: " << *its[level - 1] << " != " << cell / (size_t)pow(2, dimension) << "\n";
				++its[level - 1];
				if (its[level - 1] == refined[level - 1].end() || *its[level - 1] != cell / (size_t)pow(2, dimension)) {
					--level;
//					std::cout << "l: " << level << ", its: " << *its[level] << "\n";
					n = n >> 1;
					cell /= pow(2, dimension);
				}
			}

			if (dimension == 2) {
				D1toD2(n, cell, x, y);
			}
			if (dimension == 3) {
				D1toD3(n, cell, x, y, z);
			}

//			std::cout << "-----------------\n";
//			std::cout << "level: " << level << ", n: " << n << ", cell: " << cell << " - " << z * n * n + y * n + x << "\n";
//			std::cout << "-----\n" << n << ":" << z * n * n + y * n + x << " -> ";

			size_t xs, xe, ys, ye, zs, ze;
//			auto addNeighbors = [&] (size_t xx, size_t yy, size_t zz, size_t xs, size_t xe, size_t ys, size_t ye, size_t zs, size_t ze) {
			auto addNeighbors = [&] (size_t xx, size_t yy, size_t zz, int x, int y, int z) {
				// xx, yy -> index in deeper level
				xx = xx * (1 << (refinedxyz.size() - level - 1));
				yy = yy * (1 << (refinedxyz.size() - level - 1));
				zz = zz * (1 << (refinedxyz.size() - level - 1));
				xs = ys = zs = 0;
				xe = ye = ze = 2;
				if (x == -1) { xs = 1; }
				if (y == -1) { ys = 1; }
				if (z == -1) { zs = 1; }
				if (x == 1) { xe = 1; }
				if (y == 1) { ye = 1; }
				if (z == 1) { ze = 1; }
//				std::cout << "xx: " << xx << ", yy: " << yy << "\n";

				// l -> level, ll = divisor
				// go deeper up to my level
				size_t l = 1;
				size_t ll = N / NN;
//				std::cout << "l: " << l << ", ll: " << ll << "\n";
				while (l < level && std::binary_search(refinedxyz[l].begin(), refinedxyz[l].end(), (zz / ll) * (NN << l - 1) * (NN << l - 1) + (yy / ll) * (NN << l - 1) + xx / ll)) {
					++l;
					ll = ll >> 1;
				}

				potential.clear();
				potential.push_back(std::make_pair((NN << l - 1), (zz / ll) * (NN << l - 1) * (NN << l - 1) + (yy / ll) * (NN << l - 1) + xx / ll));
//				std::cout << "first: NN: " << NN << ", l: " << l << ", ll: " << ll << " == " << potential.back() << "\n";
				size_t nbegin = 0;
				size_t nend = 1;
				size_t xoffset = 0;
				size_t yoffset = 0;
				size_t zoffset = 0;
				// there should be more deeper levels
				while (l < refinedxyz.size()) {
					ll = ll >> 1;
					for (size_t nn = nbegin; nn < nend; ++nn) {
//						std::cout << "l: " << l << ", ll: " << ll << "\n";
						if (std::binary_search(refinedxyz[l].begin(), refinedxyz[l].end(), potential[nn].second)) {
							xoffset = 2 * (potential[nn].second % (potential[nn].first)) - xx / ll;
							yoffset = 2 * (potential[nn].second % (potential[nn].first * potential[nn].first) / (potential[nn].first)) - yy / ll;
							zoffset = 2 * (potential[nn].second / (potential[nn].first * potential[nn].first)) - zz / ll;
//							std::cout << 2 << " * " << potential[nn].second % (potential[nn].first) << " - " << xx / ll << " = " << xoffset << "\n";
//							std::cout << potential[nn] << " is IN with offset " << xoffset << ", " << yoffset << ", " << zoffset << "\n";
							for (size_t i = xs + xoffset; i < xe + xoffset; i++) {
								for (size_t j = ys + yoffset; j < ye + yoffset; j++) {
									for (size_t k = zs + zoffset; k < ze + zoffset; k++) {
//										std::cout << (zz / ll + k) << " * " << (NN << l) * (NN << l) << " = " << (zz / ll + k) * (NN << l) * (NN << l) << "\n";
//										std::cout << (yy / ll + j) << " * " << (NN << l) << " = " << (yy / ll + j) * (NN << l) << "\n";
//										std::cout << " = " << xx / ll + i << "\n";
//										std::cout << "pushback: " << std::make_pair((NN << l), (zz / ll + k) * (NN << l) * (NN << l) + (yy / ll + j) * (NN << l) + xx / ll + i) << "\n";
										potential.push_back(std::make_pair((NN << l), (zz / ll + k) * (NN << l) * (NN << l) + (yy / ll + j) * (NN << l) + xx / ll + i));
									}
								}
							}
						} else {
							neighbors.push_back(potential[nn]);
						}
					}
					nbegin = nend;
					nend = potential.size();
					++l;
				}
			};


			size_t nsize = neighbors.size();

			if (dimension == 2) {
//				D1toD2(n, cell, x, y);
				for (int ox = -1; ox <= 1; ox++) {
					for (int oy = -1; oy <= 1; oy++) {
						if (x + ox < n && y + oy < n) {
//							std::cout << x << ", " << y << ", " << z << "::" << ox << ", " << oy << ", " << 1 << "\n";
							addNeighbors(x + ox, y + oy, z, ox, oy, 1);
						}
					}
				}
			}
			if (dimension == 3) {
//				D1toD3(n, cell, x, y, z);
				for (int ox = -1; ox <= 1; ox++) {
					for (int oy = -1; oy <= 1; oy++) {
						for (int oz = -1; oz <= 1; oz++) {
							if (x + ox < n && y + oy < n && z + oz < n) {
//								std::cout << x << ", " << y << ", " << z << "::" << ox << ", " << oy << ", " << oz << "\n";
								addNeighbors(x + ox, y + oy, z + oz, ox, oy, oz);
							}
						}
					}
				}
			}


//			for (size_t i = nsize; i < neighbors.size(); i++) {
//				std::cout << neighbors[i] << " ";
//			}
//			std::cout << "\n";

			++cell;
		}
		Esutils::sortAndRemoveDuplicity(neighbors);
//		std::cout << neighbors;

		intervals.resize(neighbors.size());
		for (size_t i = 0; i < neighbors.size(); i++) {
			if (dimension == 2) {
				intervals[i].first = D2toD1(neighbors[i].first,
						neighbors[i].second % neighbors[i].first,
						neighbors[i].second / neighbors[i].first);
			}
			if (dimension == 3) {
				intervals[i].first = D3toD1(neighbors[i].first,
						neighbors[i].second % neighbors[i].first,
						neighbors[i].second % (neighbors[i].first * neighbors[i].first) / neighbors[i].first,
						neighbors[i].second / (neighbors[i].first * neighbors[i].first));
			}
			intervals[i].second = intervals[i].first + 1. / pow(neighbors[i].first, dimension);
		}

		std::sort(intervals.begin(), intervals.end());
//		std::cout << intervals;
//		for (size_t i = 0; i < intervals.size(); i++) {
//			std::cout << neighbors[i] << " -> " << intervals[i] << "\n";
//		}

		size_t unique = 0;
		for (size_t i = 1; i < intervals.size(); i++) {
			if (intervals[i].first <= intervals[unique].second + 1e-6) {
				intervals[unique].second = intervals[i].second;
			} else {
				intervals[++unique] = intervals[i];
			}
		}

		intervals.resize(unique + 1);
//		std::cout << intervals;
	});

	std::vector<int> nranks;
	for (size_t i = 0; i < intervals.size(); i++) {
		int begin = std::lower_bound(_sfcbounds.begin() + 1, _sfcbounds.end(), intervals[i].first) - _sfcbounds.begin();
		int end = std::lower_bound(_sfcbounds.begin() + 1, _sfcbounds.end(), intervals[i].second) - _sfcbounds.begin();
		for (int r = begin; r <= end; r++) {
			if (r - 1 != environment->MPIrank) {
				nranks.push_back(r - 1);
			}
		}
	}
	Esutils::sortAndRemoveDuplicity(nranks);
	Communication::serialize([&] () {
		std::cout << environment->MPIrank << "-> " << nranks;
	});


	std::vector<std::vector<eslocal> > sSize, sNodes, rSize, rNodes;
	std::vector<std::vector<EData> > sEData, rEData;
	std::vector<int> targets;

	for (int r = 0; r < environment->MPIsize; r++) {
		auto begin = std::lower_bound(permutation.begin(), permutation.end(), _sfcbounds[r], [&] (eslocal i, double b) { return partition[i] < b; });
		auto end = std::lower_bound(permutation.begin(), permutation.end(), _sfcbounds[r + 1], [&] (eslocal i, double b) { return partition[i] < b; });
		if (begin != end) {
			sSize.push_back({});
			sNodes.push_back({});
			sEData.push_back({});
			targets.push_back(r);
		}
		for (size_t n = begin - permutation.begin(); n < end - permutation.begin(); ++n) {
			sSize.back().push_back(_dMesh.esize[permutation[n]]);
			sEData.back().push_back(_dMesh.edata[permutation[n]]);
			sNodes.back().insert(sNodes.back().end(), _dMesh.enodes.begin() + edist[permutation[n]], _dMesh.enodes.begin() + edist[permutation[n] + 1]);
		}
	}

	if (!Communication::sendVariousTargets(sSize, rSize, targets)) {
		ESINFO(ERROR) << "ESPRESO internal error: distribute not sorted elements sizes.";
	}
	if (!Communication::sendVariousTargets(sEData, rEData, targets)) {
		ESINFO(ERROR) << "ESPRESO internal error: distribute not sorted element data.";
	}
	if (!Communication::sendVariousTargets(sNodes, rNodes, targets)) {
		ESINFO(ERROR) << "ESPRESO internal error: distribute not sorted element nodes.";
	}

	for (size_t r = 1; r < rSize.size(); r++) {
		rSize[0].insert(rSize[0].end(), rSize[r].begin(), rSize[r].end());
		rEData[0].insert(rEData[0].end(), rEData[r].begin(), rEData[r].end());
		rNodes[0].insert(rNodes[0].end(), rNodes[r].begin(), rNodes[r].end());
	}

	size_t enodes = 0;
	if (rSize.size()) {
		enodes = rNodes[0].size();
		permutation.resize(rSize[0].size());
		std::iota(permutation.begin(), permutation.end(), 0);
		std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) { return rEData[0][i].id < rEData[0][j].id; });

		edist = std::vector<eslocal>({ 0 });
		edist.reserve(rSize[0].size() + 1);
		for (size_t e = 0; e < rSize[0].size(); e++) {
			edist.push_back(edist.back() + rSize[0][e]);
		}
	} else {
		permutation.clear();
	}

	_dMesh.esize.clear();
	_dMesh.enodes.clear();
	_dMesh.edata.clear();
	_dMesh.esize.reserve(permutation.size());
	_dMesh.edata.reserve(permutation.size());
	_dMesh.enodes.reserve(enodes);
	for (size_t n = 0; n < permutation.size(); n++) {
		_dMesh.esize.push_back(rSize[0][permutation[n]]);
		_dMesh.edata.push_back(rEData[0][permutation[n]]);
		_dMesh.enodes.insert(_dMesh.enodes.end(), rNodes[0].begin() + edist[permutation[n]], rNodes[0].begin() + edist[permutation[n] + 1]);
	}

	_eDistribution = Communication::getDistribution(_dMesh.esize.size(), MPITools::operations().sizeToOffsetsSize_t);
}

void BalancedLoader::polishSFC()
{
	size_t threads = environment->OMP_NUM_THREADS;

	Point origin = _mesh.nodes->min, size = _mesh.nodes->max - _mesh.nodes->min + Point(1e-6, 1e-6, 1e-6);
	size_t N = 1 << (size_t)std::ceil(std::log2(environment->MPIsize));

	std::vector<eslocal> toexchange;
	std::vector<std::vector<double> > centers(threads);

//	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		centers[t].resize(_mesh.elements->dimension * (_mesh.elements->distribution[t + 1] - _mesh.elements->distribution[t]));
		Point center;
		size_t eindex = 0;
		double D1;
		size_t ncount;
		for (auto e = _mesh.elements->nodes->cbegin(t); e != _mesh.elements->nodes->cend(t); ++e, ++eindex) {
			center = Point();
			ncount = 0;
			size_t NN = N;
			while (true) {
				const Point &p = _mesh.nodes->coordinates->datatarray()[e->front()];
				size_t x = std::floor(NN * (p.x - origin.x) / size.x);
				size_t y = std::floor(NN * (p.y - origin.y) / size.y);
				D1 = D2toD1(NN, x, y);
				if (_sfcbounds[environment->MPIrank] <= D1 && D1 < _sfcbounds[environment->MPIrank + 1]) {
					break;
				}
				NN = NN << 1;
			}
			for (auto n = e->begin(); n != e->end(); ++n) {
				const Point &p = _mesh.nodes->coordinates->datatarray()[*n];
				center += p;
				size_t x = std::floor(NN * (p.x - origin.x) / size.x);
				size_t y = std::floor(NN * (p.y - origin.y) / size.y);
				D1 = D2toD1(NN, x, y);
				if (_sfcbounds[environment->MPIrank] <= D1 && D1 < _sfcbounds[environment->MPIrank + 1]) {
					++ncount;
				}
			}
			if (ncount < _mesh.elements->dimension) {
				toexchange.push_back(eindex);
			}
			center /= e->size();
			centers[t][_mesh.elements->dimension * eindex + 0] = center.x;
			centers[t][_mesh.elements->dimension * eindex + 1] = center.y;
			if (_mesh.elements->dimension == 3) {
				centers[t][_mesh.elements->dimension * eindex + 2] = center.z;
			}
		}
	}

	_mesh.elements->centers = new serializededata<eslocal, double>(_mesh.elements->dimension, centers);
}

void BalancedLoader::sortElementsVariousTargets()
{
	TimeEval time("SORT BY VARIOUS TARGETS");
	time.totalTime.startWithBarrier();

	TimeEvent e1("PREPARE");
	e1.start();

	std::vector<eslocal> edist = { 0 };
	edist.reserve(_dMesh.esize.size() + 1);
	for (size_t e = 0; e < _dMesh.esize.size(); e++) {
		edist.push_back(edist.back() + _dMesh.esize[e]);
	}

	std::vector<eslocal> permutation(_dMesh.edata.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) { return _dMesh.enodes[edist[i]] < _dMesh.enodes[edist[j]]; });

	std::vector<std::vector<eslocal> > sSize, sNodes, rSize, rNodes;
	std::vector<std::vector<EData> > sEData, rEData;
	std::vector<int> targets;

	for (int r = 0; r < environment->MPIsize; r++) {
		auto begin = std::lower_bound(permutation.begin(), permutation.end(), _nDistribution[r], [&] (eslocal i, const size_t &ID) { return _dMesh.enodes[edist[i]] < ID; });
		auto end = std::lower_bound(permutation.begin(), permutation.end(), _nDistribution[r + 1], [&] (eslocal i, const size_t &ID) { return _dMesh.enodes[edist[i]] < ID; });
		if (begin != end) {
			sSize.push_back({});
			sNodes.push_back({});
			sEData.push_back({});
			targets.push_back(r);
		}
		for (size_t n = begin - permutation.begin(); n < end - permutation.begin(); ++n) {
			sSize.back().push_back(_dMesh.esize[permutation[n]]);
			sEData.back().push_back(_dMesh.edata[permutation[n]]);
			sNodes.back().insert(sNodes.back().end(), _dMesh.enodes.begin() + edist[permutation[n]], _dMesh.enodes.begin() + edist[permutation[n] + 1]);
		}
	}

	e1.end();
	time.addEvent(e1);

	TimeEvent e2("EXCHANGE");
	e2.start();

	if (!Communication::sendVariousTargets(sSize, rSize, targets)) {
		ESINFO(ERROR) << "ESPRESO internal error: distribute not sorted elements sizes.";
	}
	if (!Communication::sendVariousTargets(sEData, rEData, targets)) {
		ESINFO(ERROR) << "ESPRESO internal error: distribute not sorted element data.";
	}
	if (!Communication::sendVariousTargets(sNodes, rNodes, targets)) {
		ESINFO(ERROR) << "ESPRESO internal error: distribute not sorted element nodes.";
	}

	for (size_t r = 1; r < rSize.size(); r++) {
		rSize[0].insert(rSize[0].end(), rSize[r].begin(), rSize[r].end());
		rEData[0].insert(rEData[0].end(), rEData[r].begin(), rEData[r].end());
		rNodes[0].insert(rNodes[0].end(), rNodes[r].begin(), rNodes[r].end());
	}

	e2.end();
	time.addEvent(e2);

	TimeEvent e3("POSTPROCESS");
	e3.start();

	size_t enodes = 0;
	size_t esize = 0;
	if (rSize.size()) {
		enodes = rNodes[0].size();
		esize = rSize[0].size();
	}

	_dMesh.esize.clear();
	_dMesh.enodes.clear();
	_dMesh.edata.clear();
	_dMesh.esize.reserve(esize);
	_dMesh.edata.reserve(esize);
	_dMesh.enodes.reserve(enodes);
	for (size_t n = 0, offset = 0; n < esize; offset += rSize[0][n++]) {
		_dMesh.esize.push_back(rSize[0][n]);
		_dMesh.edata.push_back(rEData[0][n]);
		_dMesh.enodes.insert(_dMesh.enodes.end(), rNodes[0].begin() + offset, rNodes[0].begin() + offset + rSize[0][n]);
	}

	_eDistribution = Communication::getDistribution(_dMesh.esize.size(), MPITools::operations().sizeToOffsetsSize_t);

	e3.end();
	time.addEvent(e3);

	time.totalTime.endWithBarrier();
	time.printStatsMPI();

//	DISTRIBUTION BY METIS GEOM
//	eslocal dimension = 0;
//	switch (_configuration.physics) {
//	case PHYSICS::HEAT_TRANSFER_2D:
//	case PHYSICS::STRUCTURAL_MECHANICS_2D:
//	case PHYSICS::SHALLOW_WATER_2D:
//		dimension = 2;
//		break;
//	case PHYSICS::HEAT_TRANSFER_3D:
//	case PHYSICS::STRUCTURAL_MECHANICS_3D:
//		dimension = 3;
//		break;
//	}
//
//	eslocal node;
//	std::vector<double> coordinates(_dMesh.esize.size() * dimension);
//	std::vector<eslocal> partition(_dMesh.esize.size()), frames(_eDistribution.begin(), _eDistribution.end());
//
//	for (size_t e = 0, eoffset = 0; e < _dMesh.esize.size(); eoffset += _dMesh.esize[e++]) {
//		node = _dMesh.enodes[eoffset] - _nDistribution[environment->MPIrank];
//		coordinates[dimension * e + 0] = _dMesh.coordinates[node].x;
//		coordinates[dimension * e + 1] = _dMesh.coordinates[node].y;
//		if (dimension == 3) {
//			coordinates[dimension * e + 2] = _dMesh.coordinates[node].z;
//		}
//	}
//
//	ParMETIS::call(ParMETIS::METHOD::ParMETIS_V3_PartGeom, frames.data(), NULL, NULL, dimension, coordinates.data(), 0, NULL, NULL, partition.data());
//
//	permutation.resize(partition.size());
//	std::iota(permutation.begin(), permutation.end(), 0);
//	std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) { return partition[i] < partition[j]; });
//
//
//	edist = std::vector<eslocal>({ 0 });
//	edist.reserve(_dMesh.esize.size() + 1);
//	for (size_t e = 0; e < _dMesh.esize.size(); e++) {
//		edist.push_back(edist.back() + _dMesh.esize[e]);
//	}
//
//	sSize.clear(), sNodes.clear(), rSize.clear(), rNodes.clear();
//	sEData.clear(), rEData.clear();
//	targets.clear();
//
//	eslocal prevID = -1;
//	for (auto e = permutation.begin(); e != permutation.end(); ++e) {
//		if (partition[*e] != prevID) {
//			sSize.push_back({});
//			sNodes.push_back({});
//			sEData.push_back({});
//			targets.push_back(partition[*e]);
//			prevID = partition[*e];
//		}
//		sSize.back().push_back(_dMesh.esize[*e]);
//		sEData.back().push_back(_dMesh.edata[*e]);
//		sNodes.back().insert(sNodes.back().end(), _dMesh.enodes.begin() + edist[*e], _dMesh.enodes.begin() + edist[*e + 1]);
//	}
//
//	if (!Communication::sendVariousTargets(sSize, rSize, targets)) {
//		ESINFO(ERROR) << "ESPRESO internal error: distribute not sorted elements sizes.";
//	}
//	if (!Communication::sendVariousTargets(sEData, rEData, targets)) {
//		ESINFO(ERROR) << "ESPRESO internal error: distribute not sorted element data.";
//	}
//	if (!Communication::sendVariousTargets(sNodes, rNodes, targets)) {
//		ESINFO(ERROR) << "ESPRESO internal error: distribute not sorted element nodes.";
//	}
//
//	for (size_t r = 1; r < rSize.size(); r++) {
//		rSize[0].insert(rSize[0].end(), rSize[r].begin(), rSize[r].end());
//		rEData[0].insert(rEData[0].end(), rEData[r].begin(), rEData[r].end());
//		rNodes[0].insert(rNodes[0].end(), rNodes[r].begin(), rNodes[r].end());
//	}
//
//	permutation.resize(rSize[0].size());
//	std::iota(permutation.begin(), permutation.end(), 0);
//	std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) { return rEData[0][i].id < rEData[0][j].id; });
//
//	edist = std::vector<eslocal>({ 0 });
//	edist.reserve(rSize[0].size() + 1);
//	for (size_t e = 0; e < rSize[0].size(); e++) {
//		edist.push_back(edist.back() + rSize[0][e]);
//	}
//
//	_dMesh.esize.clear();
//	_dMesh.enodes.clear();
//	_dMesh.edata.clear();
//	_dMesh.esize.reserve(permutation.size());
//	_dMesh.edata.reserve(permutation.size());
//	_dMesh.enodes.reserve(rNodes[0].size());
//	for (size_t n = 0; n < permutation.size(); n++) {
//		_dMesh.esize.push_back(rSize[0][permutation[n]]);
//		_dMesh.edata.push_back(rEData[0][permutation[n]]);
//		_dMesh.enodes.insert(_dMesh.enodes.end(), rNodes[0].begin() + edist[permutation[n]], rNodes[0].begin() + edist[permutation[n] + 1]);
//	}
//
//	_eDistribution = Communication::getDistribution(_dMesh.esize.size(), MPITools::operations().sizeToOffsetsSize_t);
}

void BalancedLoader::sortElementsAllToAll()
{

}

void BalancedLoader::sortElementsAllToAllv()
{
	TimeEval time("SORT BY ALLTOALLV");
	time.totalTime.startWithBarrier();

	TimeEvent e1("PREPARE");
	e1.start();

	std::vector<eslocal> edist = { 0 };
	edist.reserve(_dMesh.esize.size() + 1);
	for (size_t e = 0; e < _dMesh.esize.size(); e++) {
		edist.push_back(edist.back() + _dMesh.esize[e]);
	}

	std::vector<eslocal> permutation(_dMesh.edata.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) { return _dMesh.enodes[edist[i]] < _dMesh.enodes[edist[j]]; });

	std::vector<eslocal> sSize, sNodes;
	std::vector<EData> sEData;

	std::vector<int> ssize(environment->MPIsize), rsize(environment->MPIsize);
	std::vector<int> snsize(environment->MPIsize), rnsize(environment->MPIsize);

	for (int r = 0; r < environment->MPIsize; r++) {
		auto begin = std::lower_bound(permutation.begin(), permutation.end(), _nDistribution[r], [&] (eslocal i, const size_t &ID) { return _dMesh.enodes[edist[i]] < ID; });
		auto end = std::lower_bound(permutation.begin(), permutation.end(), _nDistribution[r + 1], [&] (eslocal i, const size_t &ID) { return _dMesh.enodes[edist[i]] < ID; });
		ssize[r] = end - begin;
		snsize[r] = sNodes.size();

		for (size_t n = begin - permutation.begin(); n < end - permutation.begin(); ++n) {
			sSize.push_back(_dMesh.esize[permutation[n]]);
			sEData.push_back(_dMesh.edata[permutation[n]]);
			sNodes.insert(sNodes.end(), _dMesh.enodes.begin() + edist[permutation[n]], _dMesh.enodes.begin() + edist[permutation[n] + 1]);
		}
		snsize[r] = sNodes.size() - snsize[r];
	}

	e1.end();
	time.addEvent(e1);

	TimeEvent e2("EXCHANGE SIZES");
	e2.start();

	MPI_Alltoall(ssize.data(), 1, MPI_INT, rsize.data(), 1, MPI_INT, environment->MPICommunicator);
	MPI_Alltoall(snsize.data(), 1, MPI_INT, rnsize.data(), 1, MPI_INT, environment->MPICommunicator);

	size_t rrsize = 0;
	for (int t = 0; t < environment->MPIsize; t++) {
		rrsize += rsize[t];
	}
	size_t rrnsize = 0;
	for (int t = 0; t < environment->MPIsize; t++) {
		rrnsize += rnsize[t];
	}

	_dMesh.esize.resize(rrsize);
	_dMesh.edata.resize(rrsize);
	_dMesh.enodes.resize(rrnsize);

	e2.end();
	time.addEvent(e2);

	TimeEvent e3("EXCHANGE DATA");
	e3.start();

	if (!Communication::allToAllV(sSize, _dMesh.esize, ssize, rsize)) {
		ESINFO(ERROR) << "ESPRESO internal error: distribute not sorted elements sizes.";
	}
	if (!Communication::allToAllV(sEData, _dMesh.edata, ssize, rsize)) {
		ESINFO(ERROR) << "ESPRESO internal error: distribute not sorted element data.";
	}
	if (!Communication::allToAllV(sNodes, _dMesh.enodes, snsize, rnsize)) {
		ESINFO(ERROR) << "ESPRESO internal error: distribute not sorted element nodes.";
	}

	e3.end();
	time.addEvent(e3);

	time.totalTime.endWithBarrier();
	time.printStatsMPI();

	_eDistribution = Communication::getDistribution(_dMesh.esize.size(), MPITools::operations().sizeToOffsetsSize_t);
}

void BalancedLoader::sortElementsManual()
{
	TimeEval time("MANUAL SORT");
	time.totalTime.startWithBarrier();

	TimeEvent e1("MS PREPARE");
	e1.start();

	std::vector<eslocal> edist = { 0 };
	edist.reserve(_dMesh.esize.size() + 1);
	for (size_t e = 0; e < _dMesh.esize.size(); e++) {
		edist.push_back(edist.back() + _dMesh.esize[e]);
	}

	std::vector<eslocal> permutation(_dMesh.edata.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) { return _dMesh.enodes[edist[i]] < _dMesh.enodes[edist[j]]; });

	std::vector<eslocal> sBuffer, rBuffer;

	size_t prevsize;
	for (int r = 0; r < environment->MPIsize; r++) {
		auto begin = std::lower_bound(permutation.begin(), permutation.end(), _nDistribution[r], [&] (eslocal i, const size_t &ID) { return _dMesh.enodes[edist[i]] < ID; });
		auto end = std::lower_bound(permutation.begin(), permutation.end(), _nDistribution[r + 1], [&] (eslocal i, const size_t &ID) { return _dMesh.enodes[edist[i]] < ID; });
		prevsize = sBuffer.size();
		sBuffer.push_back(0);
		sBuffer.push_back(r);
		sBuffer.push_back(end - begin);
		sBuffer.push_back(0);

		for (size_t n = begin - permutation.begin(); n < end - permutation.begin(); ++n) {
			sBuffer.push_back(_dMesh.esize[permutation[n]]);
			sBuffer.insert(sBuffer.end(), reinterpret_cast<const eslocal*>(_dMesh.edata.data() + permutation[n]), reinterpret_cast<const eslocal*>(_dMesh.edata.data() + permutation[n] + 1));
			sBuffer.insert(sBuffer.end(), _dMesh.enodes.begin() + edist[permutation[n]], _dMesh.enodes.begin() + edist[permutation[n] + 1]);
			sBuffer[prevsize + 3] += edist[permutation[n] + 1] - edist[permutation[n]];
		}
		sBuffer[prevsize] = sBuffer.size() - prevsize;
	}

	e1.end();
	time.addEvent(e1);

	TimeEvent e3("MS EXCHANGE DATA");
	e3.start();

	if (!Communication::allToAllWithDataSizeAndTarget(sBuffer, rBuffer)) {
		ESINFO(ERROR) << "ESPRESO internal error: distribute not sorted elements sizes.";
	}

	e3.end();
	time.addEvent(e3);

	_dMesh.esize.clear();
	_dMesh.edata.clear();
	_dMesh.enodes.clear();

	size_t offset = 0;
	EData edata;
	for (int r = 0; r < environment->MPIsize; r++) {
		++offset;
		size_t esize = rBuffer[++offset];
		size_t enodes = rBuffer[++offset];
		++offset;
		for (size_t e = 0; e < esize; ++e) {
			_dMesh.esize.push_back(rBuffer[offset++]);
			memcpy(&edata, rBuffer.data() + offset, sizeof(EData));
			_dMesh.edata.push_back(edata);
			offset += sizeof(EData) / sizeof(eslocal);
			_dMesh.enodes.insert(_dMesh.enodes.end(), rBuffer.begin() + offset, rBuffer.begin() + offset + _dMesh.esize.back());
			offset += _dMesh.esize.back();
		}
	}

	time.totalTime.endWithBarrier();
	time.printStatsMPI();

	_eDistribution = Communication::getDistribution(_dMesh.esize.size(), MPITools::operations().sizeToOffsetsSize_t);
}

void BalancedLoader::fillElements()
{
	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<std::vector<eslocal> > tedist(threads);
	std::vector<std::vector<eslocal> > tnodes(threads);
	std::vector<std::vector<eslocal> > eIDs(threads), rData(threads);
	std::vector<std::vector<int> > eMat(threads), eBody(threads);
	std::vector<std::vector<Element*> > epointers(threads);

	std::vector<eslocal> edist = { 0 };
	edist.reserve(_dMesh.esize.size() + 1);
	for (size_t e = 0; e < _dMesh.esize.size(); e++) {
		edist.push_back(edist.back() + _dMesh.esize[e]);
	}

	std::vector<size_t> edistribution = tarray<Point>::distribute(threads, _dMesh.esize.size());
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		if(t == 0) {
			tedist[t].insert(tedist[t].end(), edist.begin() + edistribution[t], edist.begin() + edistribution[t + 1] + 1);
		} else {
			tedist[t].insert(tedist[t].end(), edist.begin() + edistribution[t] + 1, edist.begin() + edistribution[t + 1] + 1);
		}
		tnodes[t].insert(tnodes[t].end(), _dMesh.enodes.begin() + edist[edistribution[t]], _dMesh.enodes.begin() + edist[edistribution[t + 1]]);
		eIDs[t].resize(edistribution[t + 1] - edistribution[t]);
		std::iota(eIDs[t].begin(), eIDs[t].end(), _eDistribution[environment->MPIrank] + edistribution[t]);
		epointers[t].reserve(edistribution[t + 1] - edistribution[t]);
		eBody[t].reserve(edistribution[t + 1] - edistribution[t]);

		for (size_t e = edistribution[t]; e < edistribution[t + 1]; ++e) {
			epointers[t].push_back(&_mesh._eclasses[t][_dMesh.edata[e].etype]);
			eBody[t].push_back(_dMesh.edata[e].body);
//			if (eIDs[t][e - edistribution[t]] != _dMesh.edata[e].id) {
//				ESINFO(ERROR) << "ESPRESO Workbench parser: not implemented ordering of EBLOCK elements IDs.";
//			}
		}

		if (_configuration.input == INPUT_FORMAT::WORKBENCH && _configuration.workbench.keep_material_sets) {
			eMat[t].reserve(edistribution[t + 1] - edistribution[t]);
			for (size_t e = edistribution[t]; e < edistribution[t + 1]; ++e) {
				eMat[t].push_back(_dMesh.edata[e].material);
			}
		} else {
			eMat[t].resize(edistribution[t + 1] - edistribution[t]);
		}
	}

	switch (_mesh.configuration.physics) {
	case PHYSICS::HEAT_TRANSFER_2D:
	case PHYSICS::STRUCTURAL_MECHANICS_2D:
	case PHYSICS::SHALLOW_WATER_2D:
		_mesh.elements->dimension = 2;
		break;
	case PHYSICS::HEAT_TRANSFER_3D:
	case PHYSICS::STRUCTURAL_MECHANICS_3D:
		_mesh.elements->dimension = 3;
	}
	_mesh.elements->size = _dMesh.esize.size();
	_mesh.elements->distribution = edistribution;
	_mesh.elements->IDs = new serializededata<eslocal, eslocal>(1, eIDs);
	_mesh.elements->nodes = new serializededata<eslocal, eslocal>(tedist, tnodes);
	_mesh.elements->epointers = new serializededata<eslocal, Element*>(1, epointers);
	_mesh.elements->material = new serializededata<eslocal, int>(1, eMat);
	_mesh.elements->body = new serializededata<eslocal, int>(1, eBody);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		rData[t].resize(edistribution[t + 1] - edistribution[t]);
		std::iota(rData[t].begin(), rData[t].end(), edistribution[t]);
	}
	_mesh.elementsRegions.push_back(new ElementsRegionStore("ALL_ELEMENTS"));
	_mesh.elementsRegions.back()->elements = new serializededata<eslocal, eslocal>(1, rData);
}

void BalancedLoader::fillCoordinates()
{
	size_t threads = environment->OMP_NUM_THREADS;

	if (environment->MPIsize == 1) {
		std::vector<std::vector<Point> > tcoordinates(threads);
		std::vector<std::vector<eslocal> > nIDs(threads), rData(threads);

		std::vector<size_t> cdistribution = tarray<Point>::distribute(threads, _dMesh.coordinates.size());
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			tcoordinates[t].insert(tcoordinates[t].end(), _dMesh.coordinates.begin() + cdistribution[t], _dMesh.coordinates.begin() + cdistribution[t + 1]);
			nIDs[t].insert(nIDs[t].end(), _dMesh.nIDs.begin() + cdistribution[t], _dMesh.nIDs.begin() + cdistribution[t + 1]);
		}

		_mesh.nodes->size = _dMesh.coordinates.size();
		_mesh.nodes->distribution = cdistribution;
		_mesh.nodes->IDs = new serializededata<eslocal, eslocal>(1, nIDs);
		_mesh.nodes->coordinates = new serializededata<eslocal, Point>(1, tcoordinates);
		_mesh.nodes->ranks = new serializededata<eslocal, int>(1, tarray<int>(threads, _nDistribution.back()));

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			rData[t].resize(cdistribution[t + 1] - cdistribution[t]);
			std::iota(rData[t].begin(), rData[t].end(), cdistribution[t]);
		}
		_mesh.boundaryRegions.push_back(new BoundaryRegionStore("ALL_NODES", _mesh._eclasses));
		_mesh.boundaryRegions.back()->nodes = new serializededata<eslocal, eslocal>(1, rData);

		_mesh.neighboursWithMe.push_back(environment->MPIrank);
		return;
	}

	std::vector<eslocal> nodes(_mesh.elements->nodes->datatarray().begin(), _mesh.elements->nodes->datatarray().end());
	Esutils::sortAndRemoveDuplicity(nodes);

	std::vector<std::vector<eslocal> > sBuffer;
	std::vector<int> sRanks;

	for (int t = 0; t < environment->MPIsize; t++) {
		auto begin = std::lower_bound(nodes.begin(), nodes.end(), _nDistribution[t]);
		auto end = std::lower_bound(nodes.begin(), nodes.end(), _nDistribution[t + 1]);
		if (end - begin) {
			sBuffer.push_back(std::vector<eslocal>(begin, end));
			sRanks.push_back(t);
		}
	}

	if (!Communication::sendVariousTargets(sBuffer, _rankNodeMap, sRanks, _targetRanks)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange neighbors.";
	}

	std::vector<size_t> ndistribution = tarray<Point>::distribute(threads, _dMesh.coordinates.size());
	std::vector<std::vector<std::vector<eslocal> > > backedData(threads, std::vector<std::vector<eslocal> >(_targetRanks.size()));
	std::vector<std::vector<Point> > backedCoordinates(_targetRanks.size());

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<eslocal> ranks, ranksOffset;
		std::vector<std::vector<eslocal>::const_iterator> rPointer(_targetRanks.size());
		for (size_t r = 0; r < _targetRanks.size(); r++) {
			rPointer[r] = std::lower_bound(_rankNodeMap[r].begin(), _rankNodeMap[r].end(), _nDistribution[environment->MPIrank] + ndistribution[t]);
		}
		for (size_t n = ndistribution[t]; n < ndistribution[t + 1]; ++n) {
			ranks.clear();
			ranksOffset.clear();
			for (size_t r = 0; r < _targetRanks.size(); r++) {
				if (rPointer[r] != _rankNodeMap[r].end() && *rPointer[r] == _nDistribution[environment->MPIrank] + n) {
					ranksOffset.push_back(r);
					ranks.push_back(_targetRanks[r]);
					++rPointer[r];
				}
			}
			for (size_t r = 0; r < ranks.size(); r++) {
				backedData[t][ranksOffset[r]].push_back(ranksOffset.size());
				backedData[t][ranksOffset[r]].insert(backedData[t][ranksOffset[r]].end(), ranks.begin(), ranks.end());
			}
		}
	}

	for (size_t t = 1; t < threads; t++) {
		for (size_t r = 0; r < _targetRanks.size(); r++) {
			backedData[0][r].insert(backedData[0][r].end(), backedData[t][r].begin(), backedData[t][r].end());
		}
	}

	for (size_t r = 0; r < _targetRanks.size(); r++) {
		backedCoordinates[r].reserve(_rankNodeMap[r].size());
		for (size_t n = 0; n < _rankNodeMap[r].size(); ++n) {
			backedCoordinates[r].push_back(_dMesh.coordinates[_rankNodeMap[r][n] - _nDistribution[environment->MPIrank]]);
		}
	}

	std::vector<std::vector<eslocal> > nodeRanks(sRanks.size()), allnodes(threads);
	std::vector<std::vector<Point> > coordinates(sRanks.size());

	if (!Communication::sendVariousTargets(backedData[0], nodeRanks, _targetRanks)) {
		ESINFO(ERROR) << "ESPRESO internal error: return node ranks.";
	}
	if (!Communication::sendVariousTargets(backedCoordinates, coordinates, _targetRanks)) {
		ESINFO(ERROR) << "ESPRESO internal error: return coordinates.";
	}

	size_t csize = 0;
	for (size_t i = 0; i < coordinates.size(); i++) {
		csize += coordinates[i].size();
	}

	std::vector<size_t> distribution = tarray<Point>::distribute(threads, csize);
	std::vector<std::vector<eslocal> > tIDs(threads), rankDistribution(sRanks.size());
	std::vector<std::vector<int> > rankData(sRanks.size());
	if (rankDistribution.size()) {
		rankDistribution.front().push_back(0);
	}

	#pragma omp parallel for
	for (size_t r = 0; r < sRanks.size(); r++) {
		for (size_t n = 0; n < nodeRanks[r].size(); n += nodeRanks[r][n] + 1) {
			rankData[r].insert(rankData[r].end(), nodeRanks[r].begin() + n + 1, nodeRanks[r].begin() + n + 1 + nodeRanks[r][n]);
			rankDistribution[r].push_back(rankData[r].size());
		}
	}

	Esutils::threadDistributionToFullDistribution(rankDistribution);

	for (size_t i = threads; i < sRanks.size(); i++) {
		coordinates[threads - 1].insert(coordinates[threads - 1].end(), coordinates[i].begin(), coordinates[i].end());
		rankData[threads - 1].insert(rankData[threads - 1].end(), rankData[i].begin(), rankData[i].end());
		rankDistribution[threads - 1].insert(rankDistribution[threads - 1].end(), rankDistribution[i].begin(), rankDistribution[i].end());
	}
	for (size_t i = threads; i < sRanks.size(); i++) {
		sBuffer[threads - 1].insert(sBuffer[threads - 1].end(), sBuffer[i].begin(), sBuffer[i].end());
	}
	coordinates.resize(threads);
	sBuffer.resize(threads);
	rankData.resize(threads);
	rankDistribution.resize(threads);

	serializededata<eslocal, Point>::balance(1, coordinates, &distribution);
	serializededata<eslocal, eslocal>::balance(1, sBuffer, &distribution);
	serializededata<eslocal, int>::balance(rankDistribution, rankData, &distribution);


	_mesh.nodes->size = distribution.back();
	_mesh.nodes->distribution = distribution;
	_mesh.nodes->IDs = new serializededata<eslocal, eslocal>(1, sBuffer);
	_mesh.nodes->coordinates = new serializededata<eslocal, Point>(1, coordinates);
	_mesh.nodes->ranks = new serializededata<eslocal, int>(rankDistribution, rankData);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		allnodes[t].resize(distribution[t + 1] - distribution[t]);
		std::iota(allnodes[t].begin(), allnodes[t].end(), distribution[t]);
		Esutils::sortAndRemoveDuplicity(rankData[t]);
	}

	_mesh.boundaryRegions.push_back(new BoundaryRegionStore("ALL_NODES", _mesh._eclasses));
	_mesh.boundaryRegions.back()->nodes = new serializededata<eslocal, eslocal>(1, allnodes);

	for (size_t t = 0; t < threads; t++) {
		_mesh.neighboursWithMe.insert(_mesh.neighboursWithMe.end(), rankData[t].begin(), rankData[t].end());
	}
	Esutils::sortAndRemoveDuplicity(_mesh.neighboursWithMe);

	for (size_t n = 0; n < _mesh.neighboursWithMe.size(); n++) {
		if (_mesh.neighboursWithMe[n] != environment->MPIrank) {
			_mesh.neighbours.push_back(_mesh.neighboursWithMe[n]);
		}
	}

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (auto n = _mesh.elements->nodes->begin(t)->begin(); n != _mesh.elements->nodes->end(t)->begin(); ++n) {
			*n = std::lower_bound(_mesh.nodes->IDs->datatarray().begin(), _mesh.nodes->IDs->datatarray().end(), *n) - _mesh.nodes->IDs->datatarray().begin();
		}
	}

	Communication::serialize([&] () {
		std::cout << environment->MPIrank << ": " << _mesh.neighbours;
	});
}

void BalancedLoader::checkERegions()
{
//	std::vector<MeshERegion> bregions;
//
//	for (size_t r = 0; r < _dMesh.eregions.size(); r++) {
//		if (_dMesh.eregions[r].min < _eDistribution.back() && _eDistribution.back() < _dMesh.eregions[r].max) {
//			ESINFO(ERROR) << "ESPRESO Workbench parser error: weird element region.";
//		}
//		if (_dMesh.eregions[r].min >= _eDistribution.back()) {
//			bregions.push_back(MeshERegion(std::move(_dMesh.eregions[r])));
//			_dMesh.eregions.erase(_dMesh.eregions.begin() + r--);
//		}
//	}
//
//	size_t bsize = 0;
//	std::vector<size_t> rsize = { 0 };
//	for (size_t i = 0; i < _dMesh.bregions.size(); i++) {
//		bsize += _dMesh.bregions[i].esize.size();
//		rsize.push_back(bsize);
//	}
//
//	std::vector<size_t> fdistribution = Communication::getDistribution(bsize, MPITools::operations().sizeToOffsetsSize_t);
//
//	size_t origBSize = _dMesh.bregions.size();
//
//	for (size_t r = 0; r < bregions.size(); r++) {
//		std::vector<size_t> borders;
//		for (int t = 0; t < environment->MPIsize; t++) {
//			auto begin = std::lower_bound(bregions[r].elements.begin(), bregions[r].elements.end(), fdistribution[t] + _eDistribution.back());
//			auto end = std::lower_bound(bregions[r].elements.begin(), bregions[r].elements.end(), fdistribution[t + 1] + _eDistribution.back());
//			if (begin != end) {
//				borders.push_back(*begin);
//				borders.push_back(borders.back() + end - begin);
//			}
//		}
//
//		if (!Communication::allGatherUnknownSize(borders)) {
//			ESINFO(ERROR) << "ESPRESO internal error: gather bregion borders.";
//		}
//
//		bool onlyRename = false;
//		for (size_t br = 0; br < origBSize; br++) {
//			if (_dMesh.bregions[br].min == borders.front() && _dMesh.bregions[br].max == borders.back() - 1) {
//				_dMesh.bregions[br].name = bregions[r].name;
//				onlyRename = true;
//				break;
//			}
//		}
//		if (onlyRename) {
//			continue;
//		}
//
//		std::vector<int> tRanks;
//		std::vector<std::vector<eslocal> > sBuffer, rBuffer;
//
//		for (int t = 0; t < environment->MPIsize; t++) {
//			auto begin = std::lower_bound(bregions[r].elements.begin(), bregions[r].elements.end(), fdistribution[t] + _eDistribution.back());
//			auto end = std::lower_bound(bregions[r].elements.begin(), bregions[r].elements.end(), fdistribution[t + 1] + _eDistribution.back());
//			if (begin != end) {
//				tRanks.push_back(t);
//				sBuffer.push_back(std::vector<eslocal>(begin, end));
//			}
//		}
//
//		if (!Communication::sendVariousTargets(sBuffer, rBuffer, tRanks)) {
//			ESINFO(ERROR) << "ESPRESO internal error: send boundary region indices.";
//		}
//
//		for (size_t i = 1; i < rBuffer.size(); i++) {
//			rBuffer[0].insert(rBuffer[0].end(), rBuffer[i].begin(), rBuffer[i].end());
//		}
//
//		auto cmp = [] (EData &edata, eslocal id) {
//			return edata.id < id;
//		};
//
//		_dMesh.bregions.push_back(MeshBRegion());
//		_dMesh.bregions.back().name = bregions[r].name;
//		if (rBuffer.size() && rBuffer.front().size()) {
//			for (size_t nr = 0; nr < origBSize; nr++) {
//				if (_dMesh.bregions[nr].esize.size()) {
//					auto begin = std::lower_bound(_dMesh.bregions[nr].edata.begin(), _dMesh.bregions[nr].edata.end(), rBuffer[0].front(), cmp);
//					auto end = std::lower_bound(_dMesh.bregions[nr].edata.begin(), _dMesh.bregions[nr].edata.end(), rBuffer[0].back() + 1, cmp);
//					for (size_t i = begin - _dMesh.bregions[nr].edata.begin(), nodes = 0; i < end - _dMesh.bregions[nr].edata.begin(); nodes += _dMesh.bregions[nr].esize[i++]) {
//						_dMesh.bregions.back().edata.push_back(_dMesh.bregions[nr].edata[i]);
//						_dMesh.bregions.back().enodes.insert(_dMesh.bregions.back().enodes.end(), _dMesh.bregions[nr].enodes.begin() + nodes, _dMesh.bregions[nr].enodes.begin() + nodes + _dMesh.bregions[nr].esize[i]);
//						_dMesh.bregions.back().esize.push_back(_dMesh.bregions[nr].esize[i]);
//					}
//				}
//			}
//		}
//	}
}

void BalancedLoader::addNodeRegions()
{
//	// assume sorted nodes !!
//	size_t threads = environment->OMP_NUM_THREADS;
//
//	if (environment->MPIsize == 1) {
//		for (size_t i = 0; i < _dMesh.nregions.size(); i++) {
//			_mesh.boundaryRegions.push_back(new BoundaryRegionStore(_dMesh.nregions[i].name, _mesh._eclasses));
//
//			std::vector<size_t> distribution = tarray<eslocal>::distribute(threads, _dMesh.nregions[i].nodes.size());
//			std::vector<std::vector<eslocal> > tnodes(threads);
//			#pragma omp parallel for
//			for (size_t t = 0; t < threads; t++) {
//				tnodes[t].insert(tnodes[t].end(), _dMesh.nregions[i].nodes.begin() + distribution[t], _dMesh.nregions[i].nodes.begin() + distribution[t + 1]);
//			}
//			_mesh.boundaryRegions.back()->nodes = new serializededata<eslocal, eslocal>(1, tnodes);
//
//			#pragma omp parallel for
//			for (size_t t = 0; t < threads; t++) {
//				for (auto n = _mesh.boundaryRegions.back()->nodes->begin(t)->begin(); n != _mesh.boundaryRegions.back()->nodes->end(t)->begin(); ++n) {
//					*n = std::lower_bound(_mesh.nodes->IDs->datatarray().begin(), _mesh.nodes->IDs->datatarray().end(), *n) - _mesh.nodes->IDs->datatarray().begin();
//				}
//			}
//		}
//		return;
//	}
//
//	for (size_t i = 0; i < _dMesh.nregions.size(); i++) {
//		std::sort(_dMesh.nregions[i].nodes.begin(), _dMesh.nregions[i].nodes.end());
//
//		std::vector<std::vector<eslocal> > sBuffer, rBuffer;
//		std::vector<int> sRanks, tRanks;
//
//		for (int t = 0; t < environment->MPIsize; t++) {
//			auto begin = std::lower_bound(_dMesh.nregions[i].nodes.begin(), _dMesh.nregions[i].nodes.end(), _nDistribution[t]);
//			auto end = std::lower_bound(_dMesh.nregions[i].nodes.begin(), _dMesh.nregions[i].nodes.end(), _nDistribution[t + 1]);
//			if (end - begin) {
//				sBuffer.push_back(std::vector<eslocal>(begin, end));
//				sRanks.push_back(t);
//			}
//		}
//
//		if (!Communication::sendVariousTargets(sBuffer, rBuffer, sRanks)) {
//			ESINFO(ERROR) << "ESPRESO internal error: exchange node region.";
//		}
//
//		sBuffer.clear();
//		sBuffer.resize(_targetRanks.size());
//		for (size_t r = 1; r < rBuffer.size(); r++) {
//			rBuffer[0].insert(rBuffer[0].end(), rBuffer[r].begin(), rBuffer[r].end());
//		}
//
//		if (rBuffer.size()) {
//			#pragma omp parallel for
//			for (size_t t = 0; t < _targetRanks.size(); t++) {
//				sBuffer[t].resize(rBuffer[0].size());
//				sBuffer[t].resize(std::set_intersection(_rankNodeMap[t].begin(), _rankNodeMap[t].end(), rBuffer[0].begin(), rBuffer[0].end(), sBuffer[t].begin()) - sBuffer[t].begin());
//			}
//		}
//
//		for (size_t t = 0; t < _targetRanks.size(); t++) {
//			if (sBuffer[t].size()) {
//				tRanks.push_back(t);
//			}
//		}
//		for (size_t t = 0; t < tRanks.size(); t++) {
//			sBuffer[t].swap(sBuffer[tRanks[t]]);
//			tRanks[t] = _targetRanks[tRanks[t]];
//		}
//		sBuffer.resize(tRanks.size());
//
//		rBuffer.clear();
//		if (!Communication::sendVariousTargets(sBuffer, rBuffer, tRanks)) {
//			ESINFO(ERROR) << "ESPRESO internal error: exchange node region to targets.";
//		}
//
//		for (size_t t = threads; t < rBuffer.size(); t++) {
//			rBuffer[threads - 1].insert(rBuffer[threads - 1].end(), rBuffer[t].begin(), rBuffer[t].end());
//		}
//		rBuffer.resize(threads);
//		serializededata<eslocal, eslocal>::balance(1, rBuffer);
//
//		_mesh.boundaryRegions.push_back(new BoundaryRegionStore(_dMesh.nregions[i].name, _mesh._eclasses));
//		_mesh.boundaryRegions.back()->nodes = new serializededata<eslocal, eslocal>(1, rBuffer);
//
//		#pragma omp parallel for
//		for (size_t t = 0; t < threads; t++) {
//			for (auto n = _mesh.boundaryRegions.back()->nodes->begin(t)->begin(); n != _mesh.boundaryRegions.back()->nodes->end(t)->begin(); ++n) {
//				*n = std::lower_bound(_mesh.nodes->IDs->datatarray().begin(), _mesh.nodes->IDs->datatarray().end(), *n) - _mesh.nodes->IDs->datatarray().begin();
//			}
//		}
//	}
}

void BalancedLoader::addBoundaryRegions()
{
//	size_t threads = environment->OMP_NUM_THREADS;
//
//	if (environment->MPIsize == 1) {
//		for (size_t i = 0; i < _dMesh.bregions.size(); i++) {
//			std::vector<eslocal> edist = { 0 };
//			edist.reserve(_dMesh.bregions[i].esize.size() + 1);
//			for (size_t e = 0; e < _dMesh.bregions[i].esize.size(); e++) {
//				edist.push_back(edist.back() + _dMesh.bregions[i].esize[e]);
//			}
//
//			std::vector<std::vector<eslocal> > tedist(threads), tnodes(threads);
//			std::vector<std::vector<Element*> > epointers(threads);
//			std::vector<size_t> edistribution = tarray<Point>::distribute(threads, _dMesh.bregions[i].esize.size());
//
//			tedist.front().push_back(0);
//			#pragma omp parallel for
//			for (size_t t = 0; t < threads; t++) {
//				for (size_t n = edistribution[t]; n < edistribution[t + 1]; ++n) {
//					tnodes[t].insert(tnodes[t].end(), _dMesh.bregions[i].enodes.begin() + edist[n], _dMesh.bregions[i].enodes.begin() + edist[n + 1]);
//					epointers[t].push_back(&_mesh._eclasses[t][_dMesh.bregions[i].edata[n].etype]);
//					tedist[t].push_back(tnodes[t].size());
//				}
//			}
//
//			Esutils::threadDistributionToFullDistribution(tedist);
//
//			_mesh.boundaryRegions.push_back(new BoundaryRegionStore(_dMesh.bregions[i].name, _mesh._eclasses));
//			_mesh.boundaryRegions.back()->distribution = tarray<eslocal>::distribute(threads, epointers.front().size());
//			switch (epointers.front().front()->type) {
//			case Element::TYPE::PLANE:
//				_mesh.boundaryRegions.back()->dimension = 2;
//				break;
//			case Element::TYPE::LINE:
//				_mesh.boundaryRegions.back()->dimension = 1;
//				break;
//			default:
//				ESINFO(ERROR) << "ESPRESO Workbench parser: invalid boundary region type. Have to be 3D plane or 2D line.";
//			}
//			_mesh.boundaryRegions.back()->elements = new serializededata<eslocal, eslocal>(tedist, tnodes);
//			_mesh.boundaryRegions.back()->epointers = new serializededata<eslocal, Element*>(1, epointers);
//
//			#pragma omp parallel for
//			for (size_t t = 0; t < threads; t++) {
//				for (auto n = _mesh.boundaryRegions.back()->elements->begin(t)->begin(); n != _mesh.boundaryRegions.back()->elements->end(t)->begin(); ++n) {
//					*n = std::lower_bound(_mesh.nodes->IDs->datatarray().begin(), _mesh.nodes->IDs->datatarray().end(), *n) - _mesh.nodes->IDs->datatarray().begin();
//				}
//			}
//		}
//		return;
//	}
//
//	if (_dMesh.bregions.size()) {
//		_mesh.preprocessing->linkNodesAndElements();
//	}
//	std::vector<eslocal> edistribution = _mesh.elements->gatherElementsProcDistribution();
//
//	for (size_t i = 0; i < _dMesh.bregions.size(); i++) {
//
//		std::vector<eslocal> edist = { 0 };
//		edist.reserve(_dMesh.bregions[i].esize.size() + 1);
//		for (size_t e = 0; e < _dMesh.bregions[i].esize.size(); e++) {
//			edist.push_back(edist.back() + _dMesh.bregions[i].esize[e]);
//		}
//
//		std::vector<eslocal> permutation(edist.size() - 1);
//		std::iota(permutation.begin(), permutation.end(), 0);
//		std::sort(permutation.begin(), permutation.end(), [&] (eslocal e1, eslocal e2) {
//			return _dMesh.bregions[i].enodes[edist[e1]] < _dMesh.bregions[i].enodes[edist[e2]];
//		});
//
//		std::vector<std::vector<eslocal> > sBuffer, rBuffer;
//		std::vector<int> sRanks, tRanks;
//
//		for (int t = 0; t < environment->MPIsize; t++) {
//			auto begin = std::lower_bound(permutation.begin(), permutation.end(), _nDistribution[t], [&] (eslocal e, eslocal n) { return _dMesh.bregions[i].enodes[edist[e]] < n; });
//			auto end = std::lower_bound(permutation.begin(), permutation.end(), _nDistribution[t + 1], [&] (eslocal e, eslocal n) { return _dMesh.bregions[i].enodes[edist[e]] < n; });
//			if (begin != end) {
//				sRanks.push_back(t);
//			}
//		}
//		sBuffer.resize(sRanks.size());
//
//		#pragma omp parallel for
//		for (size_t t = 0; t < sRanks.size(); t++) {
//			auto begin = std::lower_bound(permutation.begin(), permutation.end(), _nDistribution[sRanks[t]], [&] (eslocal e, eslocal n) { return _dMesh.bregions[i].enodes[edist[e]] < n; });
//			auto end = std::lower_bound(permutation.begin(), permutation.end(), _nDistribution[sRanks[t] + 1], [&] (eslocal e, eslocal n) { return _dMesh.bregions[i].enodes[edist[e]] < n; });
//			for (auto e = begin; e != end; ++e) {
//				sBuffer[t].push_back(_dMesh.bregions[i].edata[*e].etype);
//				sBuffer[t].push_back(_dMesh.bregions[i].esize[*e]);
//				for (eslocal n = 0; n < _dMesh.bregions[i].esize[*e]; ++n) {
//					sBuffer[t].push_back(_dMesh.bregions[i].enodes[edist[*e] + n]);
//				}
//			}
//		}
//
//		if (!Communication::sendVariousTargets(sBuffer, rBuffer, sRanks)) {
//			ESINFO(ERROR) << "ESPRESO internal error: exchange node region.";
//		}
//
//		sBuffer.clear();
//		sBuffer.resize(_targetRanks.size());
//
//		for (size_t r = 0; r < rBuffer.size(); r++) {
//			std::vector<eslocal> nodes;
//			for (size_t n = 0; n < rBuffer[r].size(); n += 2 + rBuffer[r][n + 1]) {
//				nodes.clear();
//				nodes.insert(nodes.end(), rBuffer[r].begin() + n + 2, rBuffer[r].begin() + n + 2 + rBuffer[r][n + 1]);
//				std::sort(nodes.begin(), nodes.end());
//				auto nbegin = std::lower_bound(nodes.begin(), nodes.end(), _nDistribution[environment->MPIrank]);
//				auto nend = std::lower_bound(nodes.begin(), nodes.end(), _nDistribution[environment->MPIrank + 1]);
//
////				#pragma omp parallel for
//				for (size_t t = 0; t < _targetRanks.size(); t++) {
//					auto it = _rankNodeMap[t].begin();
//					bool found = true;
//					for (auto current = nbegin; found && current != nend; ++current) {
//						it = std::lower_bound(it, _rankNodeMap[t].end(), *current);
//						found = it != _rankNodeMap[t].end() && *it == *current;
//					}
//					if (found) {
//						sBuffer[t].insert(sBuffer[t].end(), rBuffer[r].begin() + n, rBuffer[r].begin() + n + 2 + rBuffer[r][n + 1]);
//					}
//				}
//			}
//		}
//
//		for (size_t t = 0; t < _targetRanks.size(); t++) {
//			if (sBuffer[t].size()) {
//				tRanks.push_back(t);
//			}
//		}
//		for (size_t t = 0; t < tRanks.size(); t++) {
//			sBuffer[t].swap(sBuffer[tRanks[t]]);
//			tRanks[t] = _targetRanks[tRanks[t]];
//		}
//		sBuffer.resize(tRanks.size());
//
//		rBuffer.clear();
//		if (!Communication::sendVariousTargets(sBuffer, rBuffer, tRanks)) {
//			ESINFO(ERROR) << "ESPRESO internal error: exchange node region to targets.";
//		}
//
//		std::vector<std::vector<eslocal> > tedist(std::max((size_t)1, rBuffer.size()), { 0 }), tnodes(std::max((size_t)1, rBuffer.size()));
//		std::vector<std::vector<Element*> > epointers(std::max((size_t)1, rBuffer.size()));
//
////		#pragma omp parallel for
//		for (size_t r = 0; r < rBuffer.size(); r++) {
//			std::vector<eslocal> nodes;
//			for (size_t n = 0; n < rBuffer[r].size(); n += 2 + rBuffer[r][n + 1]) {
//				nodes.clear();
//				nodes.insert(nodes.end(), rBuffer[r].begin() + n + 2, rBuffer[r].begin() + n + 2 + rBuffer[r][n + 1]);
//				std::sort(nodes.begin(), nodes.end());
//
//				bool found = true;
//				auto it = _mesh.nodes->IDs->datatarray().begin();
//				for (auto current = nodes.begin(); found && current != nodes.end(); ++current) {
//					it = std::lower_bound(it, _mesh.nodes->IDs->datatarray().end(), *current);
//					found = it != _mesh.nodes->IDs->datatarray().end() && *it == *current;
//				}
//				if (found) {
//					tnodes[r].insert(tnodes[r].end(), rBuffer[r].begin() + n + 2, rBuffer[r].begin() + n + 2 + rBuffer[r][n + 1]);
//					tedist[r].push_back(tnodes[r].size());
//					epointers[r].push_back(&_mesh._eclasses[0][rBuffer[r][n]]);
//				}
//			}
//		}
//
////		#pragma omp parallel for
//		for (size_t r = 0; r < rBuffer.size(); r++) {
//			for (auto n = tnodes[r].begin(); n != tnodes[r].end(); ++n) {
//				*n = std::lower_bound(_mesh.nodes->IDs->datatarray().begin(), _mesh.nodes->IDs->datatarray().end(), *n) - _mesh.nodes->IDs->datatarray().begin();
//			}
//
//			auto elinks = _mesh.nodes->elements->begin();
//			std::vector<eslocal> elements;
//			eslocal size, maxsize, first = edistribution[environment->MPIrank], last = edistribution[environment->MPIrank + 1];
//			eslocal pdist = r == 0 ? 1 : 0, pnodes = 0, pclass = 0;
//			for (size_t e = 0; e < epointers[r].size(); ++e) {
//				elements.clear();
//				for (eslocal n = tedist[r][e]; n < tedist[r][e] + epointers[r][e]->coarseNodes; ++n) {
//					elements.insert(elements.end(), (elinks + tnodes[r][n])->begin(), (elinks + tnodes[r][n])->end());
//				}
//				std::sort(elements.begin(), elements.end());
//				size = 1; maxsize = 0;
//				for (size_t ei = 1; ei < elements.size(); ++ei) {
//					if (elements[ei - 1] == elements[ei]) {
//						++size;
//					} else {
//						if (maxsize < size && first <= elements[ei - 1] && elements[ei - 1] < last) {
//							maxsize = size;
//						}
//						size = 1;
//					}
//				}
//				if (maxsize < size && first <= elements.back() && elements.back() < last) {
//					maxsize = size;
//				}
//				if (maxsize == epointers[r][e]->coarseNodes) {
//					for (eslocal n = tedist[r][e]; n < tedist[r][e + 1]; ++n) {
//						tnodes[r][pnodes++] = tnodes[r][n];
//					}
//					tedist[r][pdist++] = pnodes;
//					epointers[r][pclass++] = epointers[r][e];
//				}
//			}
//			tnodes[r].resize(pnodes);
//			tedist[r].resize(pdist);
//			epointers[r].resize(pclass);
//		}
//
//		Esutils::threadDistributionToFullDistribution(tedist);
//
//		for (size_t t = threads; t < rBuffer.size(); t++) {
//			tedist[threads - 1].insert(tedist[threads - 1].end(), tedist[t].begin(), tedist[t].end());
//			tnodes[threads - 1].insert(tnodes[threads - 1].end(), tnodes[t].begin(), tnodes[t].end());
//			epointers[threads - 1].insert(epointers[threads - 1].end(), epointers[t].begin(), epointers[t].end());
//		}
//
//		tedist.resize(threads);
//		tnodes.resize(threads);
//		epointers.resize(threads);
//
//		serializededata<eslocal, eslocal>::balance(tedist, tnodes);
//		serializededata<eslocal, Element*>::balance(1, epointers);
//
//		#pragma omp parallel for
//		for (size_t t = 1; t < threads; t++) {
//			for (size_t e = 0; e < epointers[t].size(); e++) {
//				epointers[t][e] = &_mesh._eclasses[t][epointers[t][e] - _mesh._eclasses[0]];
//			}
//		}
//
//		_mesh.boundaryRegions.push_back(new BoundaryRegionStore(_dMesh.bregions[i].name, _mesh._eclasses));
//		_mesh.boundaryRegions.back()->dimension = 2;
//		if (epointers.front().size()) {
//			switch (epointers.front().front()->type) {
//			case Element::TYPE::PLANE:
//				_mesh.boundaryRegions.back()->dimension = 2;
//				break;
//			case Element::TYPE::LINE:
//				_mesh.boundaryRegions.back()->dimension = 1;
//				break;
//			default:
//				ESINFO(ERROR) << "ESPRESO Workbench parser: invalid boundary region type. Have to be 3D plane or 2D line.";
//			}
//		}
//		int dim = _mesh.boundaryRegions.back()->dimension;
//		MPI_Allreduce(&dim, &_mesh.boundaryRegions.back()->dimension, 1, MPI_INT, MPI_MIN, environment->MPICommunicator);
//
//		_mesh.boundaryRegions.back()->elements = new serializededata<eslocal, eslocal>(tedist, tnodes);
//		_mesh.boundaryRegions.back()->epointers = new serializededata<eslocal, Element*>(1, epointers);
//		_mesh.boundaryRegions.back()->distribution = _mesh.boundaryRegions.back()->epointers->datatarray().distribution();
//	}
}

void BalancedLoader::addElementRegions()
{
//	size_t threads = environment->OMP_NUM_THREADS;
//
//	if (environment->MPIsize == 1) {
//		for (size_t i = 0; i < _dMesh.eregions.size(); i++) {
//			_mesh.elementsRegions.push_back(new ElementsRegionStore(_dMesh.eregions[i].name));
//
//			std::vector<size_t> distribution = tarray<eslocal>::distribute(threads, _dMesh.eregions[i].elements.size());
//			std::vector<std::vector<eslocal> > telements(threads);
//			#pragma omp parallel for
//			for (size_t t = 0; t < threads; t++) {
//				telements[t].insert(telements[t].end(), _dMesh.eregions[i].elements.begin() + distribution[t], _dMesh.eregions[i].elements.begin() + distribution[t + 1]);
//			}
//			_mesh.elementsRegions.back()->elements = new serializededata<eslocal, eslocal>(1, telements);
//		}
//		return;
//	}
//
//	for (size_t i = 0; i < _dMesh.eregions.size(); i++) {
//		std::vector<std::vector<eslocal> > sBuffer, rBuffer;
//		std::vector<int> sRanks, tRanks;
//
//		for (int t = 0; t < environment->MPIsize; t++) {
//			auto begin = std::lower_bound(_dMesh.eregions[i].elements.begin(), _dMesh.eregions[i].elements.end(), _eDistribution[t]);
//			auto end = std::lower_bound(_dMesh.eregions[i].elements.begin(), _dMesh.eregions[i].elements.end(), _eDistribution[t + 1]);
//			if (end - begin) {
//				sBuffer.push_back(std::vector<eslocal>(begin, end));
//				sRanks.push_back(t);
//			}
//		}
//
//		if (!Communication::sendVariousTargets(sBuffer, rBuffer, sRanks)) {
//			ESINFO(ERROR) << "ESPRESO internal error: exchange node region.";
//		}
//
//		for (size_t t = threads; t < rBuffer.size(); t++) {
//			rBuffer[threads - 1].insert(rBuffer[threads - 1].end(), rBuffer[t].begin(), rBuffer[t].end());
//		}
//		rBuffer.resize(threads);
//		serializededata<eslocal, eslocal>::balance(1, rBuffer);
//
//		_mesh.elementsRegions.push_back(new ElementsRegionStore(_dMesh.eregions[i].name));
//		_mesh.elementsRegions.back()->elements = new serializededata<eslocal, eslocal>(1, rBuffer);
//
//		#pragma omp parallel for
//		for (size_t t = 0; t < threads; t++) {
//			for (auto e = _mesh.elementsRegions.back()->elements->begin(t)->begin(); e != _mesh.elementsRegions.back()->elements->end(t)->begin(); ++e) {
//				*e -= _eDistribution[environment->MPIrank];
//			}
//		}
//	}
}
