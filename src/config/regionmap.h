
#ifndef SRC_CONFIG_REGIONMAP_H_
#define SRC_CONFIG_REGIONMAP_H_

#include "configuration.h"

#include "../mesh/store/boundaryregionstore.h"

#include <map>
#include <algorithm>
#include <iostream>

namespace espreso {

struct RegionMapBase {
	std::vector<std::string> order;

	void addRegion(const std::string &name)
	{
		if (std::find(order.begin(), order.end(), name) == order.end()) {
			order.push_back(name);
		}
	}

	virtual void addIntersection(const std::string &name, std::vector<std::string> &regions) =0;

	virtual ~RegionMapBase() {}
};

template <typename TValue>
struct RegionMap: public RegionMapBase {
	std::map<std::string, TValue> regions;
	std::map<std::string, TValue> intersections;

	void addIntersection(const std::string &name, std::vector<std::string> &regions)
	{
		intersections.emplace(std::piecewise_construct, std::forward_as_tuple(name), std::forward_as_tuple(regions, this->regions));
	}
};

}



#endif /* SRC_CONFIG_REGIONMAP_H_ */
