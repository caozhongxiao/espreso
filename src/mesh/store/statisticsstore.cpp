
#include "statisticsstore.h"

#include <limits>

using namespace espreso;

Statistics::Statistics()
{
	min = std::numeric_limits<double>::max();
	max = -std::numeric_limits<double>::max();
	avg = norm = 0;
}



