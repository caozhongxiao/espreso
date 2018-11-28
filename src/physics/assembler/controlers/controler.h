
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
enum EvaluatorParameters: int;
template <typename TEBoundaries, typename TEData> class serializededata;
template <typename TEData> class tarray;
class ElementsRegionStore;

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

	virtual std::vector<double>& getSolutionStore() = 0;

	virtual void initData() = 0;
	virtual void processSolution() = 0;

	virtual void nextTime() = 0;
	virtual void parametersChanged() = 0;

	virtual void dirichletIndices(std::vector<std::vector<eslocal> > &indices) = 0;
	virtual void dirichletValues(std::vector<double> &values) = 0;

	virtual void processElements(Matrices matrices, InstanceFiller &filler) = 0;
	virtual void processBoundary(Matrices matrices, size_t rindex, InstanceFiller &filler) = 0;

	virtual ~Controler() {}
protected:
	Controler(Mesh &mesh, const Step &step);

	struct Parameter {
		serializededata<eslocal, double> *data;
		bool isConts;

		Parameter(): data(NULL), isConts(false) {}
		~Parameter();
	};

	bool setDefault(const std::map<std::string, ECFExpression> &values, double &defaultValue);
	bool setDefault(const std::map<std::string, ECFExpressionVector> &values, Point &defaultValue);

	void updateERegions(
			const std::map<std::string, ECFExpression> &settings, tarray<double> &data,
			eslocal csize, double *cbegin, double *tbegin, double time,
			EvaluatorParameters updatedParams = static_cast<EvaluatorParameters>(0));
	void updateERegions(
			const std::map<std::string, ECFExpressionVector> &settings, tarray<double> &data,
			eslocal csize, double *cbegin, double *tbegin, double time,
			EvaluatorParameters updatedParams = static_cast<EvaluatorParameters>(0));
	void updateBRegions(
			const ECFExpression &expression, Parameter &parameter, const std::vector<size_t> &distribution,
			eslocal csize, double *cbegin, double *tbegin, double time,
			EvaluatorParameters updatedParams = static_cast<EvaluatorParameters>(0));

	void averageNodeInitilization(tarray<double> &initData, std::vector<double> &averagedData);
	void nodeValuesToElements(tarray<double> &nodeData, std::vector<double> &elementData);

	Mesh &_mesh;
	const Step &_step;

	size_t _dirichletSize;
	std::vector<size_t> _nDistribution;
};

}



#endif /* SRC_PHYSICS_ASSEMBLER_CONTROLERS_CONTROLER_H_ */
