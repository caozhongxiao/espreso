
#include "dataholder.h"

#include "../basis/logging/logging.h"
#include "../solver/generic/SparseMatrix.h"

using namespace espreso;

DataHolder::DataHolder()
{
	computeKernelCallback = [] (FETI_REGULARIZATION regularization, int scSize, esint domain, bool ortogonalCluster) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: computeKernel is empty function. Fill it in assembler.";
	};

	computeKernelFromOrigKCallback = [] (FETI_REGULARIZATION regularization, int scSize, esint domain, bool ortogonalCluster) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: computeKernelFromOrigK is empty function. Fill it in assembler.";
	};

	computeKernelsFromOrigKCallback = [] (FETI_REGULARIZATION regularization, int scSize, bool ortogonalCluster) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: getKernelFromOrigK is empty function. Fill it in assembler.";
	};

	computeKernelsCallback = [] (FETI_REGULARIZATION regularization, int scSize, bool ortogonalCluster) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: computeKernels is empty function. Fill it in assembler.";
	};

	assembleB0Callback = [] (FETI_B0_TYPE type, const std::vector<SparseMatrix> &kernels) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: assembleB0 is empty function. Fill it in assembler.";
	};
}



