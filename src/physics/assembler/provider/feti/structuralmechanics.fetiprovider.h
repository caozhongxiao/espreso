
#ifndef SRC_PHYSICS_ASSEMBLER_PROVIDER_FETI_STRUCTURALMECHANICS_FETIPROVIDER_H_
#define SRC_PHYSICS_ASSEMBLER_PROVIDER_FETI_STRUCTURALMECHANICS_FETIPROVIDER_H_

#include "fetiprovider.h"

#include "../../../../basis/containers/point.h"

namespace espreso {

struct StructuralMechanicsLoadStepConfiguration;

class StructuralMechanicsFETIProvider: public FETIProvider {

public:
	StructuralMechanicsFETIProvider(StructuralMechanicsLoadStepConfiguration &configuration);

	MatrixType getMatrixType() const;
	MatrixType getMatrixType(esint domain) const;

protected:
	void prepareRegularization();

	StructuralMechanicsLoadStepConfiguration &_configuration;

	std::vector<Point> _cCenter, _cNorm;
	std::vector<double> _cr44, _cr45, _cr46, _cr55, _cr56;
	std::vector<size_t> _cNp;

	std::vector<Point> _dCenter, _dNorm;
	std::vector<double> _dr44, _dr45, _dr46, _dr55, _dr56;
	std::vector<size_t> _dNp;

};

}



#endif /* SRC_PHYSICS_ASSEMBLER_PROVIDER_FETI_STRUCTURALMECHANICS_FETIPROVIDER_H_ */
