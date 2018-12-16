
#ifndef SRC_ASSEMBLER_INSTANCE_H_
#define SRC_ASSEMBLER_INSTANCE_H_

#include <cstddef>
#include <vector>
#include <fstream>
#include <functional>

namespace espreso {

class SparseMatrix;
class Mesh;
enum class FETI_REGULARIZATION;
enum class FETI_B0_TYPE;

enum Matrices : int {
	NONE        = 0,
	K           = 1 << 0, // Stiffness matrix
	N           = 1 << 1, // Null spaces (N1, N2)
	M           = 1 << 2, // Mass matrix (only in assembler)
	R           = 1 << 3, // Residual forces (only in assembler)
	f           = 1 << 4, // Right-hand side
	B0          = 1 << 5, // Hybrid Total FETI gluing
	B1          = 1 << 6, // Total FETI gluing
	B1c         = 1 << 7, // Simple B1c
	B1duplicity = 1 << 8, // Lambdas duplicity
	primal      = 1 << 9, // Primal solution
	dual        = 1 << 10  // Dual solution (B1 * Lambdas)
};

struct DataHolder {

	DataHolder();
	~DataHolder();

	void computeKernel(FETI_REGULARIZATION regularization, int scSize, eslocal domain, bool ortogonalCluster = false) { computeKernelCallback(regularization, scSize, domain, ortogonalCluster); }
	void computeKernelFromOrigK(FETI_REGULARIZATION regularization, int scSize, eslocal domain, bool ortogonalCluster = false) { computeKernelFromOrigKCallback(regularization, scSize, domain, ortogonalCluster); }
	void computeKernelsFromOrigK(FETI_REGULARIZATION regularization, int scSize, bool ortogonalCluster = false) { computeKernelsFromOrigKCallback(regularization, scSize, ortogonalCluster); }
	void computeKernels(FETI_REGULARIZATION regularization, int scSize, bool ortogonalCluster = false) { computeKernelsCallback(regularization, scSize, ortogonalCluster); }
	void assembleB0(FETI_B0_TYPE type, const std::vector<SparseMatrix> &kernels) { assembleB0Callback(type, kernels); }

	std::vector<SparseMatrix> origK, K, origKN1, origKN2, origRegMat, N1, N2, RegMat;
	std::vector<SparseMatrix> M;
	std::vector<std::vector<double> > R, f;

	// matrices for Hybrid FETI constraints
	std::vector<SparseMatrix> B0;
	std::vector<std::vector<esglobal> > B0subdomainsMap; // TODO: not needed

	// matrices for FETI constraints
	std::vector<SparseMatrix> B1;
	std::vector<std::vector<eslocal> > B1subdomainsMap; // TODO: not needed
	std::vector<std::vector<eslocal> > B1clustersMap; // TODO: get it directly

	std::vector<std::vector<double> > B1c, LB, B1duplicity;

	std::vector<SparseMatrix> inequality;
	std::vector<std::vector<double> > inequalityC;

	// blocks types of B1
	enum CONSTRAINT {
		DIRICHLET,
		EQUALITY_CONSTRAINTS,
		INEQUALITY_CONSTRAINTS,
	};

	std::vector<size_t> block;

	std::vector<std::vector<double> > primalSolution;
	std::vector<std::vector<double> > dualSolution;

	std::function<void(FETI_REGULARIZATION regularization, size_t scSize, bool ortogonalCluster)> computeKernelsCallback;
	std::function<void(FETI_REGULARIZATION regularization, size_t scSize, bool ortogonalCluster)> computeKernelsFromOrigKCallback;
	std::function<void(FETI_REGULARIZATION regularization, size_t scSize, size_t domain, bool ortogonalCluster)> computeKernelCallback;
	std::function<void(FETI_REGULARIZATION regularization, size_t scSize, size_t domain, bool ortogonalCluster)> computeKernelFromOrigKCallback;
	std::function<void(FETI_B0_TYPE type, const std::vector<SparseMatrix> &kernels)> assembleB0Callback;
};

inline Matrices operator~(Matrices m)
{
	return static_cast<Matrices>(~static_cast<int>(m));
}

inline Matrices operator|(Matrices m1, const Matrices &m2)
{
	return static_cast<Matrices>(static_cast<int>(m1) | static_cast<int>(m2));
}

inline Matrices& operator|=(Matrices &m1, const Matrices &m2)
{
	m1 = static_cast<Matrices>(static_cast<int>(m1) | static_cast<int>(m2));
	return m1;
}

inline Matrices operator&(Matrices m1, const Matrices &m2)
{
	return static_cast<Matrices>(static_cast<int>(m1) & static_cast<int>(m2));
}

inline Matrices& operator&=(Matrices &m1, const Matrices &m2)
{
	m1 = static_cast<Matrices>(static_cast<int>(m1) & static_cast<int>(m2));
	return m1;
}

}


#endif /* SRC_ASSEMBLER_INSTANCE_H_ */
