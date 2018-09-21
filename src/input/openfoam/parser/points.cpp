
#include "../../../basis/containers/point.h"
#include "../../../basis/containers/tarray.h"
#include "../../../config/ecf/environment.h"

#include <numeric>
#include "points.h"

using namespace espreso;

bool OpenFOAMPoints::readData(std::vector<eslocal> &nIDs, std::vector<Point> &coordinates, double scaleFactor)
{
	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<size_t> tdistribution = tarray<eslocal>::distribute(threads, end - begin);

	std::vector<std::vector<Point> > points(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<Point> tpoints;

		current = begin + tdistribution[t];
		while (*current != '(') { ++current; }
		while (current < begin + tdistribution[t + 1]) {
			tpoints.push_back(Point(1, 1, 1));

			current += 1; // skip '('
			tpoints.back().x = scaleFactor * readDouble();
			tpoints.back().y = scaleFactor * readDouble();
			tpoints.back().z = scaleFactor * readDouble();
			current += 2; // skip ')\n'
		}

		points[t].swap(tpoints);
	}

	for (size_t t = 0; t < threads; t++) {
		coordinates.insert(coordinates.end(), points[t].begin(), points[t].end());
	}

	nIDs.resize(coordinates.size());
	std::iota(nIDs.begin(), nIDs.end(), 0);
	return true;
}




