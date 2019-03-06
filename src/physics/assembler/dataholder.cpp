
#include "dataholder.h"
#include "solver/generic/SparseMatrix.h"

#include "basis/utilities/debugprint.h"
#include "basis/utilities/sysutils.h"
#include "esinfo/eslog.hpp"

#include <fstream>

using namespace espreso;

DataHolder::DataHolder()
{
	computeKernelCallback = [] (FETI_REGULARIZATION regularization, int scSize, esint domain, bool ortogonalCluster) {
		eslog::globalerror("ESPRESO internal error: computeKernel is empty function. Fill it in assembler.\n");
	};

	computeKernelFromOrigKCallback = [] (FETI_REGULARIZATION regularization, int scSize, esint domain, bool ortogonalCluster) {
		eslog::globalerror("ESPRESO internal error: computeKernel is empty function. Fill it in assembler.\n");
	};

	computeKernelsFromOrigKCallback = [] (FETI_REGULARIZATION regularization, int scSize, bool ortogonalCluster) {
		eslog::globalerror("ESPRESO internal error: getKernelFromOrigK is empty function. Fill it in assembler.\n");
	};

	computeKernelsCallback = [] (FETI_REGULARIZATION regularization, int scSize, bool ortogonalCluster) {
		eslog::globalerror("ESPRESO internal error: computeKernels is empty function. Fill it in assembler.\n");
	};

	assembleB0Callback = [] (FETI_B0_TYPE type, const std::vector<SparseMatrix> &kernels) {
		eslog::globalerror("ESPRESO internal error: assembleB0 is empty function. Fill it in assembler.\n");
	};
}

template <typename TData>
static void storeData(TData &data, const std::string &dir, const std::string &name, size_t domain) {
	if (domain < data.size()) {
		std::ofstream os(utils::prepareFile(dir, name, domain));
		os << data[domain];
	}
};

void DataHolder::store(const std::string &prefix, Matrices matrices)
{
	std::string dir = utils::debugDirectory() + "/" + prefix;
	for (size_t d = 0; d < K.size(); d++) {
		if (matrices & Matrices::K) {
			storeData(K, dir, "K", d);
		}
		if (matrices & Matrices::M) {
			storeData(M, dir, "M", d);
		}
		if (matrices & Matrices::N) {
			storeData(N1, dir, "N1", d);
			storeData(N2, dir, "N2", d);
			storeData(RegMat, dir, "RegMat", d);
		}
		if (matrices & Matrices::R) {
			storeData(R, dir, "R", d);
		}
		if (matrices & Matrices::f) {
			storeData(f, dir, "f", d);
		}
		if (matrices & Matrices::Gluing) {
			storeData(B0, dir, "B0", d);
			storeData(B1, dir, "B1", d);
			storeData(B1c, dir, "B1c", d);
			storeData(B1duplicity, dir, "B1duplicity", d);
		}
		if (matrices & Matrices::Solution) {
			storeData(primalSolution, dir, "solution", d);
		}
		if (matrices & Matrices::f) {
			storeData(dualSolution, dir, "forces", d);
		}
	}
}

