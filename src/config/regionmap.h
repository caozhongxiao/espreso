
#ifndef SRC_CONFIG_REGIONMAP_H_
#define SRC_CONFIG_REGIONMAP_H_

#include "mesh/store/boundaryregionstore.h"
#include "basis/utilities/parser.h"

#include <map>

namespace espreso {

struct RegionMapBase {
	enum class RegionIntersection {
		FIRST,
		LAST,
		SUM,
		AVERAGE,
		ERROR
	};

	std::vector<std::string> order;
	RegionIntersection regions_intersection;

	void addRegion(const std::string &name)
	{
		for (size_t i = 0; i < order.size(); i++) {
			if (StringCompare::caseInsensitiveEq(order[i], name)) {
				return;
			}
		}
		order.push_back(name);
	}

	virtual void addIntersection(const std::string &name, std::vector<std::string> &regions) =0;

	RegionMapBase(): regions_intersection(RegionIntersection::LAST) {}
	virtual ~RegionMapBase() {}
};

template <typename TValue>
struct RegionMap: public RegionMapBase {
	std::map<std::string, TValue> regions;
	std::map<std::string, TValue> intersections;

	void addIntersection(const std::string &name, std::vector<std::string> &regions)
	{
		intersections.emplace(std::piecewise_construct, std::forward_as_tuple(name), std::forward_as_tuple(regions_intersection, regions, this->regions));
	}
};

}



#endif /* SRC_CONFIG_REGIONMAP_H_ */
