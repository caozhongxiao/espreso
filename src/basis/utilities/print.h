
#ifndef SRC_BASIS_UTILITIES_PRINT_H_
#define SRC_BASIS_UTILITIES_PRINT_H_

#include "basis/containers/point.h"

#include <ostream>
#include <vector>

namespace espreso {

inline std::ostream& operator<<(std::ostream& os, const Point &p) {
	os << p.x << " " << p.y << " " << p.z;
	return os;
}

template<typename T1, typename T2>
std::ostream& operator<< (std::ostream& os, const std::pair<T1, T2> &v)
{
	os << "<" << v.first << ":" << v.second << ">";
	return os;
}

template<typename T>
std::ostream& operator<< (std::ostream& os, const std::vector<T> &v)
{
	for(size_t i = 0; i < v.size(); ++i) {
		os << v[i] << " ";
	}
	os << "\n";
	return os;
}

}

#endif /* SRC_BASIS_UTILITIES_PRINT_H_ */
