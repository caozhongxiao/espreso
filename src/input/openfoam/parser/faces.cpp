
#include "faces.h"

#include "../openfoam.h"

#include "../../../basis/containers/tarray.h"
#include "../../../config/ecf/environment.h"

using namespace espreso;

bool OpenFOAMFaces::readFaces(PlainOpenFOAMData &data)
{
	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<size_t> tdistribution = tarray<eslocal>::distribute(threads, end - begin);

	std::vector<std::vector<eslocal> > fsize(threads), fnodes(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<eslocal> tsize, tnodes;

		current = begin + tdistribution[t];
		while (*(current + 1) != '(') { ++current; }
		while (current < begin + tdistribution[t + 1]) {
			tsize.push_back(readInteger());
			current += 1; // skip '('
			for (eslocal i = 0; i < tsize.back(); ++i) {
				tnodes.push_back(readInteger());
			}
			current += 2; // skip ')\n'
		}

		fsize[t].swap(tsize);
		fnodes[t].swap(tnodes);
	}

	for (size_t t = 0; t < threads; t++) {
		data.fsize.insert(data.fsize.end(), fsize[t].begin(), fsize[t].end());
		data.fnodes.insert(data.fnodes.end(), fnodes[t].begin(), fnodes[t].end());
	}

	return true;
}

bool OpenFOAMFaces::readParents(std::vector<eslocal> &data)
{
	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<size_t> tdistribution = tarray<eslocal>::distribute(threads, end - begin);

	std::vector<std::vector<eslocal> > elements(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<eslocal> telements;

		current = begin + tdistribution[t];
		while (current < begin + tdistribution[t + 1]) {
			telements.push_back(readInteger());
			current += 1; // skip '\n'
		}

		elements[t].swap(telements);
	}

	for (size_t t = 0; t < threads; t++) {
		data.insert(data.end(), elements[t].begin(), elements[t].end());
	}

	return true;
}


