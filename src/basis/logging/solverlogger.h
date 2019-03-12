
#ifndef SRC_BASIS_LOGGING_SOLVERLOGGER_H_
#define SRC_BASIS_LOGGING_SOLVERLOGGER_H_

#include <cstddef>
#include <vector>

namespace espreso {

class ProgressLogger;

class SolverLogger {

	struct Parameter {
		const char* name;

		union Data {
			int    *ivalue;
			long   *lvalue;
			size_t *svalue;
			double *dvalue;
			bool   *bvalue;
		} data;

		enum {
			INT, LONG, SIZE, DOUBLE, BOOL
		} type;
	};

public:
	void addparam(const char* name, int &value)
	{
		_parameters.push_back(Parameter{ name, Parameter::Data{ .ivalue = &value }, Parameter::INT });
	}

	void addparam(const char* name, long &value)
	{
		_parameters.push_back(Parameter{ name, Parameter::Data{ .lvalue = &value }, Parameter::LONG });
	}

	void addparam(const char* name, long unsigned int &value)
	{
		_parameters.push_back(Parameter{ name, Parameter::Data{ .svalue = &value }, Parameter::SIZE });
	}

	void addparam(const char* name, double &value)
	{
		_parameters.push_back(Parameter{ name, Parameter::Data{ .dvalue = &value }, Parameter::DOUBLE });
	}

	void addparam(const char* name, bool &value)
	{
		_parameters.push_back(Parameter{ name, Parameter::Data{ .bvalue = &value }, Parameter::BOOL });
	}

	void printheader(ProgressLogger &logger);
	void print(ProgressLogger &logger);

	SolverLogger()
	{
		_parameters.reserve(50);
	}

protected:
	std::vector<Parameter> _parameters;
};

}



#endif /* SRC_BASIS_LOGGING_SOLVERLOGGER_H_ */
