
#ifndef SRC_PHYSICS_ASSEMBLER_COMPOSER_COMPOSER_H_
#define SRC_PHYSICS_ASSEMBLER_COMPOSER_COMPOSER_H_

#include <vector>

#include "basis/matrices/matrixtype.h"

namespace espreso {

struct DataHolder;
class Controler;
enum Matrices: int;
struct NodeData;
struct SolverParameters;
class SparseMatrix;
template <typename TEBoundaries, typename TEData> class serializededata;

struct IJ {
	esint row, column;
};

inline bool operator==(const IJ &left, const IJ &right)
{
	return left.row == right.row && left.column == right.column;
}

inline bool operator!=(const IJ &left, const IJ &right)
{
	return !(left == right);
}

inline bool operator<(const IJ &left, const IJ &right)
{
	return left.row == right.row ? left.column < right.column : left.row < right.row;
}

class Composer {

public:
	Composer(Controler &controler);

	virtual void initDOFs() = 0;
	virtual void buildPatterns() = 0;
	virtual void buildDirichlet() = 0;

	virtual void initData();
	virtual void assemble(Matrices matrices, const SolverParameters &parameters) = 0;

	virtual void nextTime();
	virtual void parametersChanged();

	virtual void setDirichlet() = 0;
	virtual void synchronize() = 0;

	virtual void fillSolution() = 0;
	virtual void processSolution();

	virtual NodeData* RHS() =0;
	virtual void keepK();
	virtual void KplusAlfaM(double alfa) =0;
	virtual void applyOriginalK(NodeData *result, NodeData *x);
	virtual void applyM(NodeData *result, NodeData *x);
	virtual void enrichRHS(double alfa, NodeData* x) =0;
	virtual void RHSMinusR() =0;
	virtual void DirichletMinusRHS() =0;
	virtual void sum(NodeData *z, double alfa, NodeData* a, double beta, NodeData *b);
	virtual double multiply(NodeData *x, NodeData* y);
	virtual double residualNorm() =0;

	virtual ~Composer();

	DataHolder *data;

protected:
	virtual void apply(std::vector<SparseMatrix> &matrices, NodeData *result, NodeData *x) = 0;

	static size_t getMatrixSize(size_t size, MatrixType mtype)
	{
		switch (mtype) {
		case MatrixType::REAL_SYMMETRIC_INDEFINITE:
		case MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE:
			return (size * size - size) / 2 + size;
				break;
		case MatrixType::REAL_UNSYMMETRIC:
		default:
			return size * size;
		}
	}

	void insertKPattern(IJ *target, esint *begin, esint *end, MatrixType mtype);
	void clearMatrices(Matrices matrices, esint domain);

	Controler &_controler;

	serializededata<esint, esint> *_DOFMap;
	std::vector<esint> _dirichletMap;

	std::vector<esint> _dirichletPermutation;
};

}


#endif /* SRC_PHYSICS_ASSEMBLER_COMPOSER_COMPOSER_H_ */
