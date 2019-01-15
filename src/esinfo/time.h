
#ifndef SRC_ESINFO_TIME_H_
#define SRC_ESINFO_TIME_H_

#include <cstddef>

namespace espreso {
namespace time {
	extern size_t step;
	extern size_t substep;
	extern size_t iteration;

	extern double current;
	extern double shift; // difference between current and previous time
	extern double final;

	inline bool isInitial() { return step == 0 && substep == 0 && iteration == 0; }
	inline bool isLast() { return current == final; }
};

}


#endif /* SRC_ESINFO_TIME_H_ */
