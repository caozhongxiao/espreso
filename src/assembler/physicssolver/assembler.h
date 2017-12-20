
#ifndef SRC_ASSEMBLER_PHYSICSSOLVER_ASSEMBLER_H_
#define SRC_ASSEMBLER_PHYSICSSOLVER_ASSEMBLER_H_

#include <functional>
#include <vector>
#include <map>

namespace espreso {

struct Step;
struct Instance;
struct Physics;
class Mesh;
struct LinearSolver;
class TimeEval;
class TimeEvent;
class SparseMatrix;
class ResultStore;
enum Matrices: int;
enum class SumRestriction;

class Assembler {

public:
	Assembler(Instance &instance, Physics &physics, Mesh &mesh, Step &step, ResultStore &store, LinearSolver &linearSolver);
	~Assembler();

	void preprocessData();
	void updateMatrices(Matrices matrices);
	void processSolution();

	void setRegularizationCallback();
	void setRegularizationFromOrigKCallback();
	void setEmptyRegularizationCallback();
	void setB0Callback();

	void solve(Matrices updatedMatrices);

	void storeSolution();
	void storeSubSolution();

	void finalize();

	/// z = a * x + b + y
	void sum(std::vector<std::vector<double> > &z, double a, const std::vector<std::vector<double> > &x, double b, const std::vector<std::vector<double> > &y, const std::string &description);
	/// z = a * x + b + y (prefix variant)
	void sum(std::vector<std::vector<double> > &z, double a, const std::vector<std::vector<double> > &x, double b, const std::vector<std::vector<double> > &y, const std::vector<size_t> &prefix, const std::string &description);
	/// A += beta * B
	void sum(std::vector<SparseMatrix> &A, double beta, std::vector<SparseMatrix> &B, const std::string &description);

	/// y = A * x
	void multiply(std::vector<std::vector<double> > &y, std::vector<SparseMatrix> &A, std::vector<std::vector<double> > &x, const std::string &description);
	/// a = x * y
	double multiply(std::vector<std::vector<double> > &x, std::vector<std::vector<double> > &y, const std::string &description);

	double sumSquares(const std::vector<std::vector<double> > &data, SumRestriction restriction, const std::string &description);
	void addToDirichletInB1(double a, const std::vector<std::vector<double> > &x);
	double maxAbsValue(const std::vector<std::vector<double> > &v, const std::string &description);
	double lineSearch(const std::vector<std::vector<double> > &U, std::vector<std::vector<double> > &deltaU, std::vector<std::vector<double> > &F_ext);
	void keepK();

	Instance &instance;
	Physics &physics;
	Mesh &mesh;
	Step &step;
	ResultStore &store;
	LinearSolver &linearSolver;

protected:
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



#endif /* SRC_ASSEMBLER_PHYSICSSOLVER_ASSEMBLER_H_ */
