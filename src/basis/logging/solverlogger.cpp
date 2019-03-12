
#include "solverlogger.h"
#include "progresslogger.h"

using namespace espreso;

void SolverLogger::printheader(ProgressLogger &logger)
{
	logger.info("HEADER\n");
	for (size_t i = 0; i < _parameters.size(); i++) {
		logger.info("%s\n", _parameters[i].name);
	}
}

void SolverLogger::print(ProgressLogger &logger)
{
	for (size_t i = 0; i < _parameters.size(); i++) {
		switch (_parameters[i].type) {
		case Parameter::INT: logger.info("%d; ", *_parameters[i].data.ivalue); break;
		case Parameter::LONG: logger.info("%ld; ", *_parameters[i].data.lvalue); break;
		case Parameter::SIZE: logger.info("%ld; ", *_parameters[i].data.svalue); break;
		case Parameter::DOUBLE: logger.info("%e; ", *_parameters[i].data.dvalue); break;
		case Parameter::BOOL: logger.info("%c; ", *_parameters[i].data.bvalue ? "." : " "); break;
		}
	}
	logger.info("\n");
}
