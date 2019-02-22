
#include "timelogger.h"
#include "progressloggger.h"

#include "omp.h"
#include <cstdio>

using namespace espreso;

double TimeLogger::time()
{
	return omp_get_wtime();
}

void TimeLogger::evaluate(ProgressLogger &logger)
{
	logger.info("\n ============ TIME STATISTICS ============== \n");
	int depth = 0, width = 45;
	double prev = _init;
	for (size_t i = 0; i < _events.size(); i++) {
		switch (_events[i].type) {
		case Event::START:
			logger.info("%*s%-*s %fs\n", ++depth, " ", --width, _events[i].name, _events[i].data.time - prev);
			prev = _events[i].data.time;
			break;
		case Event::CHECKPOINT:
			logger.info("%*s%-*s %fs\n", depth, " ", width, _events[i].name, _events[i].data.time - prev);
			prev = _events[i].data.time;
			break;
		case Event::END:
			logger.info("%*s%-*s %fs\n", depth--, " ", width++, _events[i].name, _events[i].data.time - prev);
			prev = _init;
			break;
		case Event::DOUBLE:
//			logger.info("%*s%s=%f\n", depth, " ", _events[i].name, _events[i].data.dvalue);
			break;
		case Event::INT:
//			logger.info("%*s%s=%d\n", depth, " ", _events[i].name, _events[i].data.ivalue);
			break;
		case Event::LONG:
//			logger.info("%*s%s=%ld\n", depth, " ", _events[i].name, _events[i].data.lvalue);
			break;
		case Event::SIZE:
//			logger.info("%*s%s=%ld\n", depth, " ", _events[i].name, _events[i].data.svalue);
			break;
		default:
			break;
		}
	}
}



