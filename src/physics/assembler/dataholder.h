
#ifndef SRC_ASSEMBLER_INSTANCE_H_
#define SRC_ASSEMBLER_INSTANCE_H_

#include <cstddef>
#include <vector>
#include <functional>

namespace espreso {

class SparseMatrix;
class Mesh;
enum class FETI_REGULARIZATION;
enum class FETI_B0_TYPE;

enum Matrices : int {
	NONE        = 0,
	K           = 1 << 0, // Stiffness matrix
	C           = 1 << 1, // Stiffness matrix
	N           = 1 << 2, // Null spaces (N1, N2)
	M           = 1 << 3, // Mass matrix (only in assembler)
	R           = 1 << 4, // Residual forces (only in assembler)
	f           = 1 << 5, // Right-hand side
	Dirichlet   = 1 << 6, // Dirichlet boundary condition
	Gluing      = 1 << 7, // Contact boundary conditions
	Solution    = 1 << 8, // Linear System Solution
	Reactions   = 1 << 10, // Reaction forces
};

struct DataHolder {

	DataHolder();

	void store(const std::string &prefix, Matrices matrices);

	void computeKernel(FETI_REGULARIZATION regularization, int scSize, esint domain, bool ortogonalCluster = false) { computeKernelCallback(regularization, scSize, domain, ortogonalCluster); }
	void computeKernelFromOrigK(FETI_REGULARIZATION regularization, int scSize, esint domain, bool ortogonalCluster = false) { computeKernelFromOrigKCallback(regularization, scSize, domain, ortogonalCluster); }
	void computeKernelsFromOrigK(FETI_REGULARIZATION regularization, int scSize, bool ortogonalCluster = false) { computeKernelsFromOrigKCallback(regularization, scSize, ortogonalCluster); }
	void computeKernels(FETI_REGULARIZATION regularization, int scSize, bool ortogonalCluster = false) { computeKernelsCallback(regularization, scSize, ortogonalCluster); }
	void assembleB0(FETI_B0_TYPE type, const std::vector<SparseMatrix> &kernels) { assembleB0Callback(type, kernels); }

	std::vector<SparseMatrix> origK, solverK, K, origKN1, origKN2, origRegMat, N1, N2, RegMat;
	std::vector<SparseMatrix> M, C;
	std::vector<std::vector<double> > R, f, origF, solverF;

	// matrices for Hybrid FETI constraints
	std::vector<SparseMatrix> B0;
	std::vector<std::vector<esint> > B0subdomainsMap; // TODO: not needed

	// matrices for FETI constraints
	std::vector<SparseMatrix> B1;
	std::vector<std::vector<esint> > B1subdomainsMap; // TODO: not needed
	std::vector<std::vector<esint> > B1clustersMap; // TODO: get it directly

	std::vector<std::vector<double> > B1c, LB, B1duplicity;

	std::vector<SparseMatrix> inequality;
	std::vector<std::vector<double> > inequalityC;

	// blocks types of B1
	enum CONSTRAINT {
		DIRICHLET,
		EQUALITY_CONSTRAINTS,
		INEQUALITY_CONSTRAINTS,
	};

	std::vector<esint> block;

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
