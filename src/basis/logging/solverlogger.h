
#ifndef SRC_BASIS_LOGGING_SOLVERLOGGER_H_
#define SRC_BASIS_LOGGING_SOLVERLOGGER_H_

#include <cstddef>
#include <vector>

namespace espreso {

class ProgressLogger;

class SolverLogger {

	struct Parameter {
		const char* name;
		const char* shortcut;
		const char* format;
		char* csv;

		union Data {
			int         *ivalue;
			long        *lvalue;
			size_t      *svalue;
			double      *dvalue;
			const char* *cvalue;
			bool        *bvalue;
		} data;

		enum {
			INT, LONG, SIZE, DOUBLE, STRING, BOOL
		} type;
	};

public:
	void addparam(const char* name, const char* shortcut, const char* format, int &value)
	{
		_parameters.push_back(Parameter{ name, shortcut, format, NULL, Parameter::Data{ .ivalue = &value }, Parameter::INT });
	}

	void addparam(const char* name, const char* shortcut, const char* format, long &value)
	{
		_parameters.push_back(Parameter{ name, shortcut, format, NULL, Parameter::Data{ .lvalue = &value }, Parameter::LONG });
	}

	void addparam(const char* name, const char* shortcut, const char* format, long unsigned int &value)
	{
		_parameters.push_back(Parameter{ name, shortcut, format, NULL, Parameter::Data{ .svalue = &value }, Parameter::SIZE });
	}

	void addparam(const char* name, const char* shortcut, const char* format, double &value)
	{
		_parameters.push_back(Parameter{ name, shortcut, format, NULL, Parameter::Data{ .dvalue = &value }, Parameter::DOUBLE });
	}

	void addparam(const char* name, const char* shortcut, const char* format, const char* &value)
	{
		_parameters.push_back(Parameter{ name, shortcut, format, NULL, Parameter::Data{ .cvalue = &value }, Parameter::STRING });
	}

	void addparam(const char* name, const char* shortcut, bool &value)
	{
		_parameters.push_back(Parameter{ name, shortcut, NULL, NULL, Parameter::Data{ .bvalue = &value }, Parameter::BOOL });
	}

	void printheader(ProgressLogger &logger);
	void print(ProgressLogger &logger);

	SolverLogger();
	~SolverLogger();

protected:
	std::vector<Parameter> _parameters;
};

}



#endif /* SRC_BASIS_LOGGING_SOLVERLOGGER_H_ */
