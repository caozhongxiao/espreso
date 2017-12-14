
#ifndef SRC_ASSEMBLER_PHYSICS_PHYSICS_H_
#define SRC_ASSEMBLER_PHYSICS_PHYSICS_H_

#include <cstddef>
#include <string>
#include <vector>

namespace espreso {

struct Step;
enum Matrices : int;
enum class MatrixType;
template<typename TIndices> class SparseVVPMatrix;
class DenseMatrix;
class Mesh;
class Instance;
class EqualityConstraints;
class SparseMatrix;
struct PhysicsConfiguration;

template <typename TType> class edata;
struct NodeData;
struct ElementData;
struct BoundaryRegionStore;

enum class FETI_REGULARIZATION;

enum class SumOperation {
	SUM,
	AVERAGE
};

enum class SumRestriction {
	NONE,
	DIRICHLET,
	NON_DIRICHLET
};

struct Physics {
	friend class APITestESPRESODataProvider;

	Physics();
	Physics(const std::string &name, Mesh *mesh, Instance *instance, const PhysicsConfiguration *configuration);
	const std::string& name() const { return _name; }

	virtual void prepare() {};
	virtual void preprocessData(const Step &step) =0;

	virtual void updateMatrix(const Step &step, Matrices matrices);
	virtual void updateMatrix(const Step &step, Matrices matrices, size_t domain);

	virtual MatrixType getMatrixType(const Step &step, size_t domain) const =0;

	virtual void processElement(const Step &step, eslocal domain, Matrices matrices, eslocal eindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const =0;
	virtual void processFace(const Step &step, eslocal domain, const BoundaryRegionStore *region, Matrices matrices, eslocal findex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const =0;
	virtual void processEdge(const Step &step, eslocal domain, const BoundaryRegionStore *region, Matrices matrices, eslocal eindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const =0;
	virtual void processNode(const Step &step, eslocal domain, const BoundaryRegionStore *region, Matrices matrices, eslocal nindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const =0;
	virtual void processSolution(const Step &step) =0;

	virtual void makeStiffnessMatricesRegular(FETI_REGULARIZATION regularization, size_t scSize, bool ortogonalCluster);
	virtual void makeStiffnessMatrixRegular(FETI_REGULARIZATION regularization, size_t scSize, size_t domains, bool ortogonalCluster);
	virtual void analyticRegularization(size_t domain, bool ortogonalCluster) =0;

	virtual void assembleB1(const Step &step, bool withRedundantMultipliers, bool withGluing, bool withScaling);
	virtual void updateDirichletInB1(const Step &step, bool withRedundantMultipliers);
	virtual void assembleB0FromCorners();
	virtual void assembleB0FromKernels(const std::vector<SparseMatrix> &kernels);

	virtual double sumSquares(const std::vector<std::vector<double> > &data, SumOperation operation, SumRestriction restriction = SumRestriction::NONE, size_t loadStep = 0) const;

	virtual ~Physics();

protected:
	virtual void fillDOFsIndices(edata<const eslocal> &nodes, eslocal domain, std::vector<eslocal> &DOFs) const;
	virtual void insertElementToDomain(
			SparseVVPMatrix<eslocal> &K, SparseVVPMatrix<eslocal> &M,
			const std::vector<eslocal> &DOFs,
			const DenseMatrix &Ke, const DenseMatrix &Me, const DenseMatrix &Re, const DenseMatrix &fe,
			const Step &step, size_t domain, bool isBoundaryCondition);

	virtual void assembleBoundaryConditions(SparseVVPMatrix<eslocal> &K, SparseVVPMatrix<eslocal> &M, const Step &step, Matrices matrices, size_t domain);

	static void smoothstep(double &smoothStep, double &derivation, double edge0, double edge1, double value, size_t order);

	std::string _name;
	Mesh *_mesh;
	Instance *_instance;
	EqualityConstraints *_equalityConstraints;

	const PhysicsConfiguration *_configuration;
};

}

#endif /* SRC_ASSEMBLER_PHYSICS_PHYSICS_H_ */
