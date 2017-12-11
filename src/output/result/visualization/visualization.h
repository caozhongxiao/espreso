
#ifndef SRC_OUTPUT_RESULT_VISUALIZATION_VISUALIZATION_H_
#define SRC_OUTPUT_RESULT_VISUALIZATION_VISUALIZATION_H_

#include "../resultstore.h"

#include "../../../basis/containers/point.h"

namespace espreso {

struct Visualization: public ResultStoreBase {

	Visualization(const Mesh &mesh): ResultStoreBase(mesh) {}

	static Point shrink(const Point &p, const Point &ccenter, const Point &dcenter, double cratio, double dratio) {
		Point point = ccenter + (p - ccenter) * cratio;
		point = dcenter + (point - dcenter) * dratio;
		return point;
	}
};

}


#endif /* SRC_OUTPUT_RESULT_VISUALIZATION_VISUALIZATION_H_ */
