
#ifndef SRC_PHYSICS_ASSEMBLER_COMPOSER_COMPOSER_H_
#define SRC_PHYSICS_ASSEMBLER_COMPOSER_COMPOSER_H_

#include <vector>

#include "../../../basis/matrices/matrixtype.h"

namespace espreso {

class Controler;
enum Matrices: int;
struct NodeData;
struct SolverParameters;
template <typename TEBoundaries, typename TEData> class serializededata;

struct IJ {
	eslocal row, column;
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

	virtual void initData();
	virtual void initDirichlet() = 0;
	virtual void assemble(Matrices matrices, const SolverParameters &parameters) = 0;

	virtual void nextTime();
	virtual void parametersChanged();

	virtual void setDirichlet() = 0;
	virtual void synchronize() = 0;

	virtual void fillSolution() = 0;
	virtual void processSolution();

	NodeData* RHS();

	void keepK();
	virtual void KplusAlfaM(double alfa) =0;
	virtual void applyM(NodeData *y, NodeData *x) =0;
	virtual void applyOriginalK(NodeData *y, NodeData *x) =0;
	virtual void enrichRHS(double alfa, NodeData* a) =0;
	virtual void RHSMinusR() =0;
	virtual void DirichletMinusRHS() =0;
	void sum(NodeData *z, double alfa, NodeData* a, double beta, NodeData *b);
	double multiply(NodeData *x, NodeData* y);
	virtual double residualNorm() =0;

	virtual ~Composer() {}

protected:
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

	void insertKPattern(IJ *target, eslocal *begin, eslocal *end, MatrixType mtype);
	void clearMatrices(Matrices matrices, eslocal domain);

	Controler &_controler;

	serializededata<eslocal, eslocal> *_DOFMap;
	std::vector<eslocal> _dirichletMap;

	std::vector<eslocal> _dirichletPermutation;
};

}


#endif /* SRC_PHYSICS_ASSEMBLER_COMPOSER_COMPOSER_H_ */
