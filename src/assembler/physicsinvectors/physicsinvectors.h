
#ifndef SRC_ASSEMBLER_PHYSICSINVECTORS_PHYSICSINVECTORS_H_
#define SRC_ASSEMBLER_PHYSICSINVECTORS_PHYSICSINVECTORS_H_

#include <cstddef>
#include <string>
#include <vector>

namespace espreso {

struct Step;
class Mesh;
class Instance;
class DenseMatrix;
struct PhysicsConfiguration;
template <typename TEBoundaries, typename TEData> class serializededata;

struct DOFs {
	int element, faces, edges, nodes, midnodes;
};

struct PhysicsInVectors {

	PhysicsInVectors(const std::string &name, Mesh &mesh, Instance &instance, Step &step, const PhysicsConfiguration &configuration);
	const std::string& name() const { return _name; }

	void buildCSRPattern();
	void buildGlobalCSRPattern();

	void computeValues();
	void computeGlobalValues();
	void synchronize();

	virtual void initData() =0;
	virtual void updateData() = 0;

	virtual void setDirichlet() =0;

	virtual ~PhysicsInVectors() {};

	// deprecated functions
	virtual eslocal processElement(eslocal eindex, eslocal nindex, DenseMatrix &Ke, DenseMatrix &fe) =0;

protected:
	std::string _name;
	Mesh &_mesh;
	Instance &_instance;
	Step &_step;

	const PhysicsConfiguration &_configuration;

	DOFs _DOFs;

	std::vector<size_t> _nDistribution;
	serializededata<eslocal, eslocal>* _domainNodes;
	std::vector<std::vector<eslocal> > _DOFsPermutation;

	size_t _localK, _localRHS;
	std::vector<eslocal> _pK, _pRHS, _neighRHSSize, _globalIndices;
};

}


#endif /* SRC_ASSEMBLER_PHYSICSINVECTORS_PHYSICSINVECTORS_H_ */
