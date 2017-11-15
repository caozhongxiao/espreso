
#ifndef SRC_ASSEMBLER_PHYSICS_PRECOMPUTED_H_
#define SRC_ASSEMBLER_PHYSICS_PRECOMPUTED_H_

#include "physics.h"

namespace espreso {

struct Precomputed: public virtual Physics
{
	enum SolutionIndex {
		DECOMPOSED = 0,
		MERGED     = 1,

		SIZE       = 2
	};

	Precomputed(Mesh *mesh, Instance *instance, MatrixType type, double *rhs, size_t rhsSize);

	std::vector<size_t> solutionsIndicesToStore() const;
	std::vector<std::pair<ElementType, Property> > propertiesToStore() const;

	MatrixType getMatrixType(const Step &step, size_t domain) const;
	bool isMatrixTimeDependent(const Step &step) const;
	bool isMatrixTemperatureDependent(const Step &step) const;

	void preprocessData(const Step &step);

	void analyticRegularization(size_t domain, bool ortogonalCluster);

	void updateMatrix(const Step &step, Matrices matrices, size_t domain, const std::vector<Solution*> &solution);
	void assembleB0FromCorners();

	void processElement(const Step &step, Matrices matrices, eslocal eindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const;
	void processFace(const Step &step, Matrices matrices, eslocal findex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const;
	void processEdge(const Step &step, Matrices matrices, eslocal eindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const;
	void processNode(const Step &step, Matrices matrices, eslocal nindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const;
	void processSolution(const Step &step);

	virtual ~Precomputed() {}

protected:
	MatrixType _mtype;
	double *_rhs;
	size_t _rhsSize;
	std::vector<std::vector<double> > _mergedSolution;
};

}



#endif /* SRC_ASSEMBLER_PHYSICS_PRECOMPUTED_H_ */
