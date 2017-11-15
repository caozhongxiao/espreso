
#ifndef SRC_ASSEMBLER_PHYSICS_STRUCTURALMECHANICS_H_
#define SRC_ASSEMBLER_PHYSICS_STRUCTURALMECHANICS_H_

#include "physics.h"
#include "../../basis/point/point.h"

namespace espreso {

struct StructuralMechanicsConfiguration;
struct ResultsSelectionConfiguration;

struct StructuralMechanics: public virtual Physics
{
	StructuralMechanics(const StructuralMechanicsConfiguration &configuration, const ResultsSelectionConfiguration &propertiesConfiguration);

	virtual std::vector<size_t> solutionsIndicesToStore() const;

	virtual MatrixType getMatrixType(const Step &step, size_t domain) const;
	virtual bool isMatrixTimeDependent(const Step &step) const;
	virtual bool isMatrixTemperatureDependent(const Step &step) const;
	virtual void prepare();
	virtual void preprocessData(const Step &step);

protected:
	enum SolutionIndex: size_t {
		DISPLACEMENT = 0,

		SIZE         = 1
	};

	static size_t offset;

	const StructuralMechanicsConfiguration &_configuration;
	const ResultsSelectionConfiguration &_propertiesConfiguration;

	// to handle with non-continuous partition
	std::vector<Point> _cCenter, _cNorm;
	std::vector<double> _cr44, _cr45, _cr46, _cr55, _cr56;
	std::vector<size_t> _cNp;

	std::vector<Point> _dCenter, _dNorm;
	std::vector<double> _dr44, _dr45, _dr46, _dr55, _dr56;
	std::vector<size_t> _dNp;
};

}


#endif /* SRC_ASSEMBLER_PHYSICS_STRUCTURALMECHANICS_H_ */
