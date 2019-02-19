
#include "timelogger.h"

#include "omp.h"
#include <cstdio>

using namespace espreso;

double TimeLogger::time()
{
	return omp_get_wtime();
}

void TimeLogger::evaluate()
{
	for (size_t i = 0; i < _events.size(); i++) {
		switch (_events[i].type) {
		case Event::START:
		case Event::CHECKPOINT:
		case Event::END:
			printf("%s=%f\n", _events[i].name, _events[i].data.time);
			break;
		default:
			break;
		}
	}
}



