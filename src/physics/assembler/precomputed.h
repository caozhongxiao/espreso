
#ifndef SRC_ASSEMBLER_PHYSICS_PRECOMPUTED_H_
#define SRC_ASSEMBLER_PHYSICS_PRECOMPUTED_H_

#include "../assembler/physics.h"

namespace espreso {

struct Precomputed: public virtual Physics
{
	enum SolutionIndex {
		DECOMPOSED = 0,
		MERGED     = 1,

		SIZE       = 2
	};

	Precomputed(Mesh *mesh, Instance *instance, MatrixType type, double *rhs, size_t rhsSize);

	MatrixType getMatrixType(size_t domain) const;

	void preprocessData();
	void setDirichlet();

	void analyticRegularization(size_t domain, bool ortogonalCluster);

	void updateMatrix(Matrices matrices, size_t domain);
	void assembleB0FromCorners();

	void processBEM(eslocal domain, Matrices matrices);
	void processElement(eslocal domain, Matrices matrices, eslocal eindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const;
	void processFace(eslocal domain, const BoundaryRegionStore *region, Matrices matrices, eslocal findex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const;
	void processEdge(eslocal domain, const BoundaryRegionStore *region, Matrices matrices, eslocal eindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const;
	void processNode(eslocal domain, const BoundaryRegionStore *region, Matrices matrices, eslocal nindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const;
	void processSolution();

	virtual ~Precomputed() {}

protected:
	MatrixType _mtype;
	double *_rhs;
	size_t _rhsSize;
	std::vector<std::vector<double> > _mergedSolution;
};

}



#endif /* SRC_ASSEMBLER_PHYSICS_PRECOMPUTED_H_ */
