
#include "../mesh/mesh.h"
#include "../mesh/store/elementstore.h"
#include "../solver/generic/SparseMatrix.h"
#include "../basis/logging/logging.h"
#include "dataholder.h"

using namespace espreso;

DataHolder::DataHolder()
{
//	domainDOFCount.resize(domains);
//
//	origK.resize(domains);
//	origKN1.resize(domains);
//	origKN2.resize(domains);
//	origRegMat.resize(domains);
//
//	K.resize(domains);
//	N1.resize(domains);
//	N2.resize(domains);
//	RegMat.resize(domains);
//
//	M.resize(domains);
//	R.resize(domains);
//	f.resize(domains);
//
//	B0.resize(domains);
//	B0subdomainsMap.resize(domains);
//
//	B1.resize(domains);
//	B1subdomainsMap.resize(domains);
//	B1duplicity.resize(domains);
//	B1c.resize(domains);
//	LB.resize(domains);
//
//	inequality.resize(domains);
//	inequalityC.resize(domains);

//	block.resize(3);
//
//	for (size_t d = 0; d < domains; d++) {
//		B0[d].rows = 0;
//		B0[d].cols = 0;
//		B0[d].nnz = 0;
//		B0[d].type = 'G';
//
//		B1[d].rows = 0;
//		B1[d].cols = 0;
//		B1[d].nnz = 0;
//		B1[d].type = 'G';
//	}

	computeKernelCallback = [] (FETI_REGULARIZATION regularization, size_t scSize, size_t domain, bool ortogonalCluster) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: computeKernel is empty function. Fill it in assembler.";
	};

	computeKernelFromOrigKCallback = [] (FETI_REGULARIZATION regularization, size_t scSize, size_t domain, bool ortogonalCluster) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: computeKernelFromOrigK is empty function. Fill it in assembler.";
	};

	computeKernelsFromOrigKCallback = [] (FETI_REGULARIZATION regularization, size_t scSize, bool ortogonalCluster) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: getKernelFromOrigK is empty function. Fill it in assembler.";
	};

	computeKernelsCallback = [] (FETI_REGULARIZATION regularization, size_t scSize, bool ortogonalCluster) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: computeKernels is empty function. Fill it in assembler.";
	};

	assembleB0Callback = [] (FETI_B0_TYPE type, const std::vector<SparseMatrix> &kernels) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: assembleB0 is empty function. Fill it in assembler.";
	};
}

DataHolder::~DataHolder()
{

}



