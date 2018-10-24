
#include "mpi.h"
#include "spacefillingcurve.h"

#include "../../basis/containers/tarray.h"
#include "../../basis/logging/logging.h"
#include "../../basis/utilities/utils.h"
#include "../../config/ecf/environment.h"

#include <utility>
#include <limits>
#include <algorithm>

using namespace espreso;

SpaceFillingCurve::SpaceFillingCurve(size_t dimension, size_t depth, std::vector<Point> &coordinates)
: _dimension(dimension), _depth(depth), _n(1 << depth), _n2(_n * _n), _n3(_n * _n * _n), _refinedsfc(1, std::vector<size_t>(1))
{
	if (_dimension != 2 && _dimension != 3) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: incorrect mesh dimension ='" << _dimension << "'.";
	}
	size_t threads = environment->OMP_NUM_THREADS;

	double dmax = std::numeric_limits<double>::max();
	std::vector<Point> mins(threads, Point(dmax, dmax, dmax)), maxs(threads, Point(-dmax, -dmax, -dmax));

	std::vector<size_t> cdistribution = tarray<size_t>::distribute(threads, coordinates.size());

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		Point tmin = Point(dmax, dmax, dmax), tmax = Point(-dmax, -dmax, -dmax);
		for (size_t i = cdistribution[t]; i < cdistribution[t + 1]; i++) {
			tmin.x = std::min(tmin.x, coordinates[i].x);
			tmin.y = std::min(tmin.y, coordinates[i].y);
			tmin.z = std::min(tmin.z, coordinates[i].z);
			tmax.x = std::max(tmax.x, coordinates[i].x);
			tmax.y = std::max(tmax.y, coordinates[i].y);
			tmax.z = std::max(tmax.z, coordinates[i].z);
		}

		mins[t] = tmin;
		maxs[t] = tmax;
	}

	for (size_t t = 1; t < threads; t++) {
		mins[0].x = std::min(mins[t].x, mins[0].x);
		mins[0].y = std::min(mins[t].y, mins[0].y);
		mins[0].z = std::min(mins[t].z, mins[0].z);
		maxs[0].x = std::max(maxs[t].x, maxs[0].x);
		maxs[0].y = std::max(maxs[t].y, maxs[0].y);
		maxs[0].z = std::max(maxs[t].z, maxs[0].z);
	}

	MPI_Allreduce(&mins[0].x, &_origin.x, 3, MPI_DOUBLE, MPI_MIN, environment->MPICommunicator);
	MPI_Allreduce(&maxs[0].x, &_size.x, 3, MPI_DOUBLE, MPI_MAX, environment->MPICommunicator);

	_size -= _origin - Point(1e-6, 1e-6, 1e-6);
}

void SpaceFillingCurve::finishLevel(size_t depth)
{
	Esutils::sortAndRemoveDuplicity(_refinedsfc[depth]);
}

void SpaceFillingCurve::SCFToXYZ()
{
	size_t x, y, z = 0, n = 1;
	for (size_t i = 0; i < _refinedsfc.size(); ++i, n = n << 1) {
		_refinedxyz.push_back({});
		for (size_t j = 0; j < _refinedsfc[i].size(); j++) {
			if (_dimension == 2) {
				D1toD2(n, _refinedsfc[i][j], x, y);
				_refinedxyz[i].push_back(y * n + x);
			}
			if (_dimension == 3) {
				D1toD3(n, _refinedsfc[i][j], x, y, z);
				_refinedxyz[i].push_back(z * n * n + y * n + x);
			}
		}
		std::sort(_refinedxyz[i].begin(), _refinedxyz[i].end());
	}
}

size_t SpaceFillingCurve::buckets(size_t depth) const
{
	return std::pow(1 << depth, _dimension);
}

void SpaceFillingCurve::iterateBuckets(size_t begin, size_t end, std::function<void(size_t depth, size_t index)> callback) const
{
	std::vector<std::vector<size_t>::const_iterator> rindices;

	size_t bsize = bucketSize();
	size_t coarsenig = buckets(depth());

	for (size_t i = 0, c = coarsenig; i < _refinedsfc.size(); i++, c /= bsize) {
		rindices.push_back(std::lower_bound(_refinedsfc[i].begin(), _refinedsfc[i].end(), begin / c));
	}

	size_t depth, step;
	for (size_t bucket = begin; bucket < end; bucket += step) {
		for (size_t i = 0, c = coarsenig; i < _refinedsfc.size(); i++, c /= bsize) {
			if (rindices[i] != _refinedsfc[i].end() && *rindices[i] < bucket / c) {
				++rindices[i];
			}
		}
		depth = 0;
		step = coarsenig;
		while (depth < _depth && rindices[depth] != _refinedsfc[depth].end() && *rindices[depth] == bucket / step) {
			++depth;
			step /= bsize;
		}

		callback(depth, bucket / step);
	}
}

void SpaceFillingCurve::addSFCNeighbors(size_t depth, size_t index, std::vector<std::pair<size_t, size_t> > &neighbors)
{
	size_t x, y, z = 0, nsize;
	if (_dimension == 2) {
		D1toD2(1 << depth, index, x, y);
		nsize = neighbors.size();
		addXYNeighbors(depth, x, y, neighbors);
		for (size_t i = nsize; i < neighbors.size(); i++) {
			size_t n = 1 << neighbors[i].first;
			neighbors[i].second = D2toD1(n, neighbors[i].second % n, neighbors[i].second / n);
		}
	}
	if (_dimension == 3) {
		D1toD3(1 << depth, index, x, y, z);
		nsize = neighbors.size();
		addXYZNeighbors(depth, x, y, z, neighbors);
		for (size_t i = nsize; i < neighbors.size(); i++) {
			size_t n = 1 << neighbors[i].first;
			neighbors[i].second = D3toD1(n, neighbors[i].second % n, neighbors[i].second % (n * n) / n, neighbors[i].second / (n * n));
		}
	}
}

std::pair<size_t, size_t> SpaceFillingCurve::getXYZBucket(size_t depth, size_t x, size_t y, size_t z)
{
	std::pair<size_t, size_t> bucket;

	size_t coarsenig = 1 << depth;
	size_t cx, cy, cz;

	size_t d = 0, n = 1;
	do {
		cx = x / coarsenig;
		cy = y / coarsenig;
		cz = z / coarsenig;

		bucket.first = d;
		bucket.second = cz * n * n + cy * n + cx;

		coarsenig /= 2; n *= 2;
	} while (++d <= depth && std::binary_search(_refinedxyz[bucket.first].begin(), _refinedxyz[bucket.first].end(), bucket.second));

	return bucket;
}

void SpaceFillingCurve::addXYNeighbors(size_t depth, size_t x, size_t y, std::vector<std::pair<size_t, size_t> > &neighbors)
{
	std::vector<std::pair<size_t, size_t> > potential;

	size_t n = 1 << depth;
	size_t nn, xx, yy;

	for (int ox = -1; ox <= 1; ox++) {
		for (int oy = -1; oy <= 1; oy++) {
			if (x + ox < n && y + oy < n) {

				potential.clear();
				potential.push_back(getXYZBucket(depth, x + ox, y + oy, 0));

				// neighbors should be finer
				for (size_t i = 0; i < potential.size(); i++) {
					if (potential[i].first < _refinedxyz.size() && std::binary_search(_refinedxyz[potential[i].first].begin(), _refinedxyz[potential[i].first].end(), potential[i].second)) {
						nn = 1 << potential[i].first;
						xx = potential[i].second % nn;
						yy = potential[i].second / nn;
						nn = nn << 1;
						for (size_t oox = (ox == -1 ? 1 : 0); oox < (ox == 1 ? 1 : 2); oox++) {
							for (size_t ooy = (oy == -1 ? 1 : 0); ooy < (oy == 1 ? 1 : 2); ooy++) {
								potential.push_back(std::pair<size_t, size_t>(potential[i].first + 1, nn * (2 * yy + ooy) + 2 * xx + oox));
							}
						}
					} else {
						neighbors.push_back(potential[i]);
					}
				}

			}
		}
	}
}

void SpaceFillingCurve::addXYZNeighbors(size_t depth, size_t x, size_t y, size_t z, std::vector<std::pair<size_t, size_t> > &neighbors)
{
	std::vector<std::pair<size_t, size_t> > potential;

	size_t n = 1 << depth;
	size_t nn, xx, yy, zz;

	for (int ox = -1; ox <= 1; ox++) {
		for (int oy = -1; oy <= 1; oy++) {
			for (int oz = -1; oz <= 1; oz++) {
				if (x + ox < n && y + oy < n && z + oz < n) {

					potential.clear();
					potential.push_back(getXYZBucket(depth, x + ox, y + oy, z + oz));

					for (size_t i = 0; i < potential.size(); i++) {
						if (potential[i].first + 1 < _refinedxyz.size() && std::binary_search(_refinedxyz[potential[i].first].begin(), _refinedxyz[potential[i].first].end(), potential[i].second)) {
							nn = 1 << potential[i].first;
							xx = potential[i].second % nn;
							yy = potential[i].second % (nn * nn) / nn;
							zz = potential[i].second / (nn * nn);
							nn = nn << 1;
							for (size_t oox = (ox == -1 ? 1 : 0); oox < (ox == 1 ? 1 : 2); oox++) {
								for (size_t ooy = (oy == -1 ? 1 : 0); ooy < (oy == 1 ? 1 : 2); ooy++) {
									for (size_t ooz = (oz == -1 ? 1 : 0); ooz < (oz == 1 ? 1 : 2); ooz++) {
										potential.push_back(std::pair<size_t, size_t>(potential[i].first + 1, nn * nn * (2 * zz + ooz) + nn * (2 * yy + ooy) + 2 * xx + oox));
									}
								}
							}
						} else {
							neighbors.push_back(potential[i]);
						}
					}

				}
			}
		}
	}
}




