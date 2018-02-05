
#ifndef SRC_OUTPUT_RESULT_VISUALIZATION_DISTRIBUTEDVISUALIZATION_H_
#define SRC_OUTPUT_RESULT_VISUALIZATION_DISTRIBUTEDVISUALIZATION_H_

#include "../resultstore.h"

#include "../../../basis/containers/point.h"

namespace espreso {

struct DistributedVisualization: public ResultStoreBase {

	DistributedVisualization(const Mesh &mesh): ResultStoreBase(mesh) {}

	virtual bool isCollected() { return false; }
	virtual bool isDistributed() { return true; }

	static Point shrink(const Point &p, const Point &ccenter, const Point &dcenter, double cratio, double dratio) {
		Point point = ccenter + (p - ccenter) * cratio;
		point = dcenter + (point - dcenter) * dratio;
		return point;
	}
};

}


#endif /* SRC_OUTPUT_RESULT_VISUALIZATION_DISTRIBUTEDVISUALIZATION_H_ */
