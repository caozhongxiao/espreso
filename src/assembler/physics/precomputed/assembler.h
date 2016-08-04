
#ifndef SRC_ASSEMBLER_PHYSICS_PRECOMPUTED_ASSEMBLER_H_
#define SRC_ASSEMBLER_PHYSICS_PRECOMPUTED_ASSEMBLER_H_

#include "../assembler.h"

namespace espreso {

struct PrecomputedPhysics: public Physics {

	virtual bool singular() const
	{
		return true;
	}

	virtual void assemble()
	{
		ESINFO(PROGRESS2) << "Assemble matrices K and RHS.";
		K.resize(_mesh.parts());
		f.resize(_mesh.parts());
		cilk_for (size_t p = 0; p < _mesh.parts(); p++) {
			composeSubdomain(p);
			K[p].mtype = mtype;
			ESINFO(PROGRESS2) << Info::plain() << ".";
		}
		ESINFO(PROGRESS2);
	}

	PrecomputedPhysics(
			APIMesh &mesh, SparseMatrix::MatrixType mtype,
			const std::vector<Property> elementDOFs,
			const std::vector<Property> faceDOFs,
			const std::vector<Property> edgeDOFs,
			const std::vector<Property> pointDOFs,
			const std::vector<Property> midPointDOFs,
			double *rhs, size_t rhs_size)
	: Physics(mesh, mtype, elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs), _apimesh(mesh), rhs(rhs), rhs_size(rhs_size) {};
	virtual ~PrecomputedPhysics() {};

protected:
	virtual void composeSubdomain(size_t subdomain) =0;

	const APIMesh &_apimesh;
	double *rhs;
	size_t rhs_size;
};

}


#endif /* SRC_ASSEMBLER_PHYSICS_PRECOMPUTED_ASSEMBLER_H_ */
