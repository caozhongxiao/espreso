
#ifndef SRC_PHYSICS_PROVIDER_COLLECTIVEPROVIDER_H_
#define SRC_PHYSICS_PROVIDER_COLLECTIVEPROVIDER_H_

#include "provider.h"

namespace espreso {

class CollectiveProvider: public Provider {

public:
	CollectiveProvider(Instance &instance, Composer &composer, Mesh &mesh, Step &step, ResultStore &store, LinearSolver &linearSolver);
	~CollectiveProvider();

	void setRegularizationCallback();
	void setRegularizationFromOrigKCallback();
	void setEmptyRegularizationCallback();
	void setB0Callback();

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
};

}


#endif /* SRC_PHYSICS_PROVIDER_COLLECTIVEPROVIDER_H_ */
