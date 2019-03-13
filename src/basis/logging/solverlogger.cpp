
#include "solverlogger.h"
#include "progresslogger.h"

#include <cstring>

using namespace espreso;

SolverLogger::SolverLogger()
{

}

SolverLogger::~SolverLogger()
{
	for (size_t i = 0; i < _parameters.size(); i++) {
		if (_parameters[i].csv) {
			delete[] _parameters[i].csv;
		}
	}
}

void SolverLogger::printheader(ProgressLogger &logger)
{
	const char* prefix = " ";
	const char* suffix = "; ";

	logger.info(" ## ESPRESO PHYSICAL SOLVER\n");
	for (size_t i = 0; i < _parameters.size(); i++) {
		logger.info(" ## %s\n", _parameters[i].name);
	}
	logger.info(" ##########################\n  ");
	for (size_t i = 0; i < _parameters.size(); i++) {
		logger.info(" %s; ", _parameters[i].shortcut);

		_parameters[i].csv = new char[strlen(_parameters[i].format) + strlen(prefix) + strlen(suffix)];
		strcpy(_parameters[i].csv, prefix);
		strcat(_parameters[i].csv, _parameters[i].format);
		strcat(_parameters[i].csv, suffix);
	}
	logger.info("\n");
}

void SolverLogger::print(ProgressLogger &logger)
{
	logger.info("  ");
	for (size_t i = 0; i < _parameters.size(); i++) {
		switch (_parameters[i].type) {
		case Parameter::INT: logger.info(_parameters[i].csv, *_parameters[i].data.ivalue); break;
		case Parameter::LONG: logger.info(_parameters[i].csv, *_parameters[i].data.lvalue); break;
		case Parameter::SIZE: logger.info(_parameters[i].csv, *_parameters[i].data.svalue); break;
		case Parameter::DOUBLE: logger.info(_parameters[i].csv, *_parameters[i].data.dvalue); break;
		case Parameter::STRING: logger.info(_parameters[i].csv, *_parameters[i].data.cvalue); break;
		case Parameter::BOOL: logger.info(_parameters[i].csv, *_parameters[i].data.bvalue ? "." : " "); break;
		}
	}
	logger.info("\n");
}
