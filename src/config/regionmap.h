
#ifndef SRC_CONFIG_REGIONMAP_H_
#define SRC_CONFIG_REGIONMAP_H_

#include "configuration.h"

#include <map>
#include <algorithm>

namespace espreso {

template <typename TValue>
struct RegionMap {
	std::map<std::string, TValue> regions;

	std::vector<std::string> order;
	bool areDisjunkt;

	void addRegion(const std::string &name)
	{
		if (std::find(order.begin(), order.end(), name) == order.end()) {
			order.push_back(name);
		}
	}

	RegionMap(): areDisjunkt(true) {}
};

}



#endif /* SRC_CONFIG_REGIONMAP_H_ */
