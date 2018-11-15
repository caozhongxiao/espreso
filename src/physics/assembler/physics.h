
#ifndef SRC_ASSEMBLER_PHYSICS_PHYSICS_H_
#define SRC_ASSEMBLER_PHYSICS_PHYSICS_H_

#include <cstddef>
#include <string>
#include <vector>

namespace bem4i { template< class LO, class SC > struct bem4iData; }

namespace espreso {

struct Step;
enum Matrices : int;
enum class MatrixType;
template<typename TIndices> class SparseVVPMatrix;
class DenseMatrix;
class Mesh;
class Instance;
class Constraints;
class SparseMatrix;
struct PhysicsConfiguration;

template <typename TType> class edata;
template <typename TEBoundaries, typename TEData> class serializededata;
struct NodeData;
struct ElementData;
struct BoundaryRegionStore;

enum class FETI_REGULARIZATION;

enum class SumRestriction {
	NONE,
	DIRICHLET,
	NON_DIRICHLET
};

struct Physics {
	friend class APITestESPRESODataProvider;

	Physics();
	Physics(const std::string &name, Mesh *mesh, Instance *instance, Step *step, const PhysicsConfiguration *configuration, int DOFs);
	const std::string& name() const { return _name; }


//public:
	virtual void initLocalDOFs(std::vector<eslocal> &offsets) {}
	virtual void initGlobalDOFs(eslocal &offset) {}

	virtual void buildLocalCSRPattern() {}
	virtual void buildGlobalCSRPattern() {}

	virtual void initData() {}

//protected:

	virtual void initLocalNodeUniformDOFs(std::vector<eslocal> &offsets, eslocal multiplier);
	virtual void initGlobalNodeUniformDOFs(eslocal &offset, eslocal multiplier);

	virtual void buildLocalNodeUniformCSRPattern(eslocal multiplier);
	virtual void buildGlobalNodeUniformCSRPattern(eslocal multiplier);

	serializededata<eslocal, eslocal> *_DOF;
	std::vector<std::vector<eslocal> > _KPermutation, _RHSPermutation;

	// OLD METHODS

	virtual void prepare() {};
	virtual void preprocessData() =0;
	virtual void setDirichlet() =0;

	virtual void updateMatrix(Matrices matrices);
	virtual void updateMatrix(Matrices matrices, size_t domain);

	virtual MatrixType getMatrixType(size_t domain) const =0;

	virtual void processBEM(eslocal domain, Matrices matrices) =0;
	virtual void processElement(eslocal domain, Matrices matrices, eslocal eindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const =0;
	virtual void processFace(eslocal domain, const BoundaryRegionStore *region, Matrices matrices, eslocal findex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const =0;
	virtual void processEdge(eslocal domain, const BoundaryRegionStore *region, Matrices matrices, eslocal eindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const =0;
	virtual void processNode(eslocal domain, const BoundaryRegionStore *region, Matrices matrices, eslocal nindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const =0;
	virtual void processSolution() =0;

	virtual void makeStiffnessMatricesRegular(FETI_REGULARIZATION regularization, size_t scSize, bool ortogonalCluster);
	virtual void makeStiffnessMatrixRegular(FETI_REGULARIZATION regularization, size_t scSize, size_t domains, bool ortogonalCluster);
	virtual void analyticRegularization(size_t domain, bool ortogonalCluster) =0;

	virtual void assembleB1(bool withRedundantMultipliers, bool withGluing, bool withScaling);
	virtual void updateDirichletInB1(bool withRedundantMultipliers);
	virtual void updateDuplicity();
	virtual void assembleB0FromCorners();
	virtual void assembleB0FromKernels(const std::vector<SparseMatrix> &kernels);

	virtual double sumSquares(const std::vector<std::vector<double> > &data, SumRestriction restriction = SumRestriction::NONE) const;

	virtual ~Physics();

//protected:

	void printInvalidElement(eslocal eindex) const;

	static void smoothstep(double &smoothStep, double &derivation, double edge0, double edge1, double value, size_t order);

	std::string _name;
	Mesh *_mesh;
	Instance *_instance;
	Step *_step;
	Constraints *_constraints;

	const PhysicsConfiguration *_configuration;

	int _DOFs;

	bool _hasBEM;
	std::vector<int> _BEMDomain;
	std::vector<bem4i::bem4iData<eslocal, double>* > _BEMData;

	mutable size_t _invalidElements;
};

}

#endif /* SRC_ASSEMBLER_PHYSICS_PHYSICS_H_ */
