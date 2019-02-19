
#include "dataholder.h"
#include "esinfo/eslog.hpp"
#include "solver/generic/SparseMatrix.h"

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



