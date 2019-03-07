
#ifndef SRC_PHYSICS_ASSEMBLER_COMPOSER_COMPOSER_H_
#define SRC_PHYSICS_ASSEMBLER_COMPOSER_COMPOSER_H_

#include <cstddef>
#include <vector>

#include "basis/matrices/matrixtype.h"

namespace espreso {

struct DataHolder;
class Controller;
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
	Composer(Controller &controler);

	virtual void initDOFs() = 0;
	virtual void buildPatterns() = 0;
	virtual void buildDirichlet() = 0;
	virtual void buildMVData() = 0;

	virtual void assemble(Matrices matrices, const SolverParameters &parameters) = 0;
	virtual void setDirichlet(Matrices matrices, double reduction, const std::vector<double> &subtraction) = 0;

	virtual void fillSolution() = 0;
	virtual void computeReactionForces() = 0;
	virtual void processSolution();

	virtual void keepK();
	virtual void keepSolverK();
	virtual void keepRHS();
	virtual void keepSolverRHS();
	virtual void KplusAlfaM(double alfa) =0;
	virtual void alfaKplusBetaM(double alfa, double beta) =0;
	virtual void applyOriginalK(NodeData* result, NodeData* x);
	virtual void applyM(NodeData* result, NodeData* x);
	virtual void enrichRHS(double alfa, NodeData* x) =0;
	virtual void enrichSolution(double alfa, NodeData* x);
	virtual void RHSMinusR();
	virtual void sum(NodeData *z, double alfa, NodeData* a, double beta, NodeData *b);
	virtual double lineSearch(NodeData *U, const SolverParameters &parameters);
	virtual double mult(NodeData *x, NodeData* y);
	virtual double residualNormNumerator() =0;
	virtual double residualNormDenominator() =0;

	virtual ~Composer();

	DataHolder *data;

protected:
	virtual void apply(std::vector<SparseMatrix> &matrices, std::vector<double> &result, std::vector<double> &x) = 0;
	virtual void gather(std::vector<double> &out, std::vector<std::vector<double> > &in) = 0;

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

	Controller &_controler;

	esint _foreignDOFs;
	serializededata<esint, esint> *_DOFMap;
	std::vector<esint> _dirichletMap;

	std::vector<esint> _dirichletPermutation;
};

}


#endif /* SRC_PHYSICS_ASSEMBLER_COMPOSER_COMPOSER_H_ */
