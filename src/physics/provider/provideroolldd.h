
#ifndef SRC_PHYSICS_PROVIDER_PROVIDEROOLLDD_H_
#define SRC_PHYSICS_PROVIDER_PROVIDEROOLLDD_H_

#include <functional>
#include <vector>
#include <map>

namespace espreso {

struct DataHolder;
struct Composer;
class Mesh;
struct LinearSolver;
class TimeEval;
class TimeEvent;
class SparseMatrix;
enum Matrices: int;
enum class SumRestriction;

class GeneralHeatTransfer2D;

class ProviderOOLLDD {

public:
	ProviderOOLLDD(DataHolder &instance, Composer &composer, Mesh &mesh, LinearSolver &linearSolver);
	virtual ~ProviderOOLLDD();

	virtual void preprocessData();
	virtual void updateStructuralMatrices(Matrices matrices);
	virtual void updateGluingMatrices(Matrices matrices);
	virtual void processSolution();

	virtual void setRegularizationCallback() {}
	virtual void setRegularizationFromOrigKCallback() {}
	virtual void setEmptyRegularizationCallback() {}
	virtual void setB0Callback() {}

	virtual void solve(Matrices updatedMatrices);

	virtual void nextTime();
	virtual void parametersChanged();

	virtual void storeSolution();
	virtual void storeSubSolution() {}

	void finalize();

	/// z = a * x + b + y
	virtual void sum(std::vector<std::vector<double> > &z, double a, const std::vector<std::vector<double> > &x, double b, const std::vector<std::vector<double> > &y, const std::string &description) {}
	/// z = a * x + b + y (prefix variant)
	virtual void sum(std::vector<std::vector<double> > &z, double a, const std::vector<std::vector<double> > &x, double b, const std::vector<std::vector<double> > &y, const std::vector<size_t> &prefix, const std::string &description) {}
	/// A += beta * B
	virtual void sum(std::vector<SparseMatrix> &A, double beta, std::vector<SparseMatrix> &B, const std::string &description) {}

	/// y = A * x
	virtual void multiply(std::vector<std::vector<double> > &y, std::vector<SparseMatrix> &A, std::vector<std::vector<double> > &x, const std::string &description) {}
	/// a = x * y
	virtual double multiply(std::vector<std::vector<double> > &x, std::vector<std::vector<double> > &y, const std::string &description) { return 0; }

	virtual double sumSquares(const std::vector<std::vector<double> > &data, SumRestriction restriction, const std::string &description) { return 0; }
	virtual void addToDirichletInB1(double a, const std::vector<std::vector<double> > &x) {}
	virtual double maxAbsValue(const std::vector<std::vector<double> > &v, const std::string &description) { return 0; }
	virtual double lineSearch(const std::vector<std::vector<double> > &U, std::vector<std::vector<double> > &deltaU, std::vector<std::vector<double> > &F_ext) { return 0; }
	virtual void keepK() {}

	DataHolder &instance;
	Composer &composer;
	Mesh &mesh;
	LinearSolver &linearSolver;

protected:
	std::string mNames(espreso::Matrices matrices, const std::string &prefix = "");
	void timeWrapper(const std::string &action, std::function<void(void)> operations);

	bool checkForStore(const std::string &name);
	void storeMatrices(Matrices matrices, size_t domain);
	void storeWrapper(const std::string &name, Matrices matrices);
	void storeWrapper(const std::string &name, Matrices matrices, size_t domain);
	void storeWrapper(const std::string &name, std::vector<SparseMatrix> &matrices);
	void storeWrapper(const std::string &name, std::vector<std::vector<double> > &data);

	TimeEval *_timeStatistics;
	std::map<std::string, TimeEvent*> _timeEvents;
};

}


#endif /* SRC_PHYSICS_PROVIDER_PROVIDEROOLLDD_H_ */
