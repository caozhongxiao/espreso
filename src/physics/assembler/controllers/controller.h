
#ifndef SRC_PHYSICS_ASSEMBLER_CONTROLLERS_CONTROLLER_H_
#define SRC_PHYSICS_ASSEMBLER_CONTROLLERS_CONTROLLER_H_

#include <cstddef>
#include <functional>
#include <map>

#include "basis/matrices/denseMatrix.h"

namespace espreso {

struct PhysicsConfiguration;
struct SolverParameters;
struct Point;
class ECFExpression;
class ECFExpressionVector;
enum Matrices: int;
struct BoundaryRegionStore;
struct NodeData;
enum EvaluatorParameters: int;
template <typename TEBoundaries, typename TEData> class serializededata;
template <typename TEData> class tarray;

class Controller
{

public:
	struct InstanceFiller {
		esint begin, end;

		DenseMatrix Ke, Me, Re, fe;

		std::function<void(size_t)> insert;
	};

	virtual NodeData* solution() = 0;
	virtual const PhysicsConfiguration& configuration() const =0;

	virtual void processSolution() = 0;

	virtual void nextTime() = 0;
	virtual void parametersChanged() = 0;

	virtual void dirichletIndices(std::vector<std::vector<esint> > &indices) = 0;
	virtual void dirichletValues(std::vector<double> &values) = 0;

	virtual void processBEMdomain(esint domain, double *values);
	virtual void fillBEMInterior(esint domain, double *values);

	virtual void processElements(Matrices matrices, const SolverParameters &parameters, InstanceFiller &filler) = 0;
	virtual void processBoundary(Matrices matrices, const SolverParameters &parameters, size_t rindex, InstanceFiller &filler) = 0;

	virtual ~Controller() {}
protected:
	Controller(int dimension);

	struct Parameter {
		serializededata<esint, double> *data;
		int dimension;
		bool isConts;

		Parameter(int dimension): data(NULL), dimension(dimension), isConts(false) {}
		~Parameter();
	};

	void setCoordinates(Parameter &coordinates, serializededata<esint, esint> *procNodes = NULL);
	void setDirichlet(std::vector<double> &solution);

	void takeKernelParam(Parameter &parameter, Parameter &previous);


	void initKernelParam(Parameter &parameter, const std::map<std::string, ECFExpressionVector> &values, double defaultValue);
	void initKernelParam(Parameter &parameter, const std::map<std::string, ECFExpression> &values, double defaultValue);
	void initKernelParam(Parameter &parameter, const std::map<std::string, ECFExpression> &values, double defaultValue,
			BoundaryRegionStore *region);
	void initKernelParam(Parameter &parameter, const ECFExpression &value, double defaultValue,
			BoundaryRegionStore *region);

	void updateKernelParam(Parameter &parameter, const std::map<std::string, ECFExpression> &values, double *cbegin, double *tbegin);
	void updateKernelParam(Parameter &parameter, const std::map<std::string, ECFExpressionVector> &values, double *cbegin, double *tbegin);
	void updateKernelParam(Parameter &parameter, const ECFExpression &value, double *cbegin, double *tbegin,
			BoundaryRegionStore *region);

	void kernelToBoundary(Parameter &parameter, Parameter &boundary, BoundaryRegionStore *region);
	void kernelToNodes(Parameter &kvalues, std::vector<double> &nvalues);
	void kernelToElements(Parameter &kvalues, std::vector<double> &evalues);
	void nodesToKernels(std::vector<double> &nvalues, Parameter &kvalues, serializededata<esint, esint> *procNodes = NULL);

	int _dimension;
	size_t _dirichletSize;
	std::vector<size_t> _nDistribution;
};

}



#endif /* SRC_PHYSICS_ASSEMBLER_CONTROLLERS_CONTROLLER_H_ */
