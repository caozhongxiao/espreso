
#ifndef SRC_GLOBALS_TIME_H_
#define SRC_GLOBALS_TIME_H_

#include <cstddef>

namespace espreso {

struct time {
	static bool isInitial() { return step == 0 && substep == 0 && iteration == 0; }
	static bool isLast() { return current == final; }

	static size_t step;
	static size_t substep;
	static size_t iteration;

	static double current;
	static double shift; // difference between current and previous time
	static double final;
};

}


#endif /* SRC_GLOBALS_TIME_H_ */
