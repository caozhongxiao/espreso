
#ifndef SRC_OUTPUT_RESULT_VISUALIZATION_SEPARATED_SEPARATEDVISUALIZATION_H_
#define SRC_OUTPUT_RESULT_VISUALIZATION_SEPARATED_SEPARATEDVISUALIZATION_H_

#include "output/result/visualization/visualization.h"

#include "basis/containers/point.h"

namespace espreso {

struct SeparatedVisualization: public Visualization {

	SeparatedVisualization(const Mesh &mesh, const OutputConfiguration &configuration): Visualization(mesh, configuration) {}

	virtual bool isCollected() { return false; }
	virtual bool isSeparated() { return true; }

	static Point shrink(const Point &p, const Point &ccenter, const Point &dcenter, double cratio, double dratio) {
		Point point = ccenter + (p - ccenter) * cratio;
		point = dcenter + (point - dcenter) * dratio;
		return point;
	}
};

}


#endif /* SRC_OUTPUT_RESULT_VISUALIZATION_SEPARATED_SEPARATEDVISUALIZATION_H_ */
