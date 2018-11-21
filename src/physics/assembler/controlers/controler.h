
#ifndef SRC_PHYSICS_ASSEMBLER_CONTROLERS_CONTROLER_H_
#define SRC_PHYSICS_ASSEMBLER_CONTROLERS_CONTROLER_H_

#include <cstddef>
#include <functional>
#include <map>

#include "../../../basis/matrices/denseMatrix.h"

namespace espreso {

enum class MatrixType;
struct Point;
class Mesh;
struct Step;
class ECFExpression;
class ECFExpressionVector;
enum Matrices: int;
template <typename TEBoundaries, typename TEData> class serializededata;
template <typename TEData> class tarray;

class Controler
{

public:
	struct InstanceFiller {
		eslocal begin, end;

		DenseMatrix Ke, Me, Re, fe;

		std::function<void(size_t)> insert;
	};

	virtual MatrixType getMatrixType() const = 0;
	virtual MatrixType getMatrixType(size_t domain) const = 0;

	virtual void initData() = 0;
	virtual void updateData() = 0;

	virtual void processElements(Matrices matrices, InstanceFiller &filler) = 0;
	virtual void processBoundary(Matrices matrices, InstanceFiller &filler) = 0;

	virtual ~Controler() {}
protected:
	Controler(Mesh &mesh, const Step &step);

	bool tryElementConstness(const std::map<std::string, ECFExpression> &values, double &defaultValue);
	bool tryElementConstness(const std::map<std::string, ECFExpressionVector> &values, Point &defaultValue);

	void evaluate(
			const std::map<std::string, ECFExpression> &settings, tarray<double> &data,
			eslocal csize, double *cbegin);
	void evaluate(
			const std::map<std::string, ECFExpressionVector> &settings, tarray<double> &data,
			eslocal csize, double *cbegin);

	void nodeValuesToElements(tarray<double> &nodeData, std::vector<double> &elementData);

	struct Parameter {
		serializededata<eslocal, double> *data;
		bool isConts;

		Parameter(): data(NULL), isConts(false) {}
		~Parameter();
	};

	Mesh &_mesh;
	const Step &_step;

	std::vector<size_t> _nDistribution;
};

}



#endif /* SRC_PHYSICS_ASSEMBLER_CONTROLERS_CONTROLER_H_ */
