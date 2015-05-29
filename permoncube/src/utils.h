
#ifndef PM_UTILS_H_
#define PM_UTILS_H_

#include "settings.h"

namespace permoncube {

template <class TElement>
class Utils {

public:
	static void globalNodesCount(const Settings &settings, size_t nodes[]);
	static void clusterNodesCount(const Settings &settings, size_t nodes[]);
	static size_t clusterElementsCount(const Settings &settings);
};

}


#include "utils.hpp"


#endif /* PM_UTILS_H_ */
