
#ifndef SRC_CONFIG_ECF_MATERIAL_COORDINATESYSTEM_H_
#define SRC_CONFIG_ECF_MATERIAL_COORDINATESYSTEM_H_

#include <config/holders/expression.h>
#include "config/configuration.h"

namespace espreso {

struct CoordinateSystemConfiguration: public ECFObject {

	enum class TYPE {
		CARTESIAN,
		CYLINDRICAL,
		SPHERICAL
	};

	TYPE type;
	DIMENSION dimension;

	ECFExpressionVector rotation;
	ECFExpressionVector center;

	CoordinateSystemConfiguration();

	void createScalingMatrix(std::vector<double> &m, double x, double y, double z=0) const;

	void createTranslationMatrixToCenter(std::vector<double> &m) const;
	void createTranslationMatrixToZero(std::vector<double> &m) const;
	void createTranslationMatrix(std::vector<double> &m, double x, double y, double z=0) const;
	void createRotationMatrix(std::vector<double> &m) const;

	void multiplyTransformationMatrices(std::vector<double> &left, std::vector<double> &result) const;
	Point applyTransformation(std::vector<double> &m, const Point &p) const;

};


}

#endif /* SRC_CONFIG_ECF_MATERIAL_COORDINATESYSTEM_H_ */
