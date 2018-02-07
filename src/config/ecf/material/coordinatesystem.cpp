
#include "coordinatesystem.h"

#include "../../configuration.hpp"

#include "../../../basis/containers/point.h"
#include "../../../basis/evaluator/evaluator.h"
#include "../../../wrappers/math/math.h"

#include "mkl.h"

espreso::CoordinateSystemConfiguration::CoordinateSystemConfiguration()
: dimension(DIMENSION::D3), rotation(dimension, ECFMetaData::getcoordinatevariables(), "0"), center(dimension, ECFMetaData::getcoordinatevariables(), "0")
{
	type = TYPE::CARTESIAN;

	REGISTER(type, ECFMetaData()
			.setdescription({ "Type" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("CARTESIAN").setdescription("Cartesian system."))
			.addoption(ECFOption().setname("CYLINDRICAL").setdescription("Cylindrical system."))
			.addoption(ECFOption().setname("SPHERICAL").setdescription("Spherical system.").allowonly([&] () { return dimension == DIMENSION::D3; })));

	REGISTER(rotation, ECFMetaData()
			.setdescription({ "A rotation of the material." })
			.allowonly([&] () { return type == TYPE::CARTESIAN; }));

	REGISTER(center, ECFMetaData()
			.setdescription({ "A center of the material." })
			.allowonly([&] () { return type != TYPE::CARTESIAN; }));
}


void espreso::CoordinateSystemConfiguration::createTranslationMatrix(std::vector<double> &m, double x, double y, double z) const {
	switch (dimension) {
			case DIMENSION::D3: {
				m.resize(16,0);
				m[0 *4 + 0] = 1;
				m[1 *4 + 1] = 1;
				m[2 *4 + 2] = 1;
				m[3 *4 + 3] = 1;
				m[0 *4 + 3] = x;
				m[1 *4 + 3] = y;
				m[2 *4 + 3] = z;

			} break;

			case DIMENSION::D2: {
				m.resize(9,0);
				m[0 *3 + 0] = 1;
				m[1 *3 + 1] = 1;
				m[2 *3 + 2] = 1;
				m[0 *3 + 2] = x;
				m[1 *3 + 2] = y;

			} break;
			default:
				ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: unsupported operation.";
		}
}
void espreso::CoordinateSystemConfiguration::createScalingMatrix(std::vector<double> &m, double x, double y, double z) const {
	switch (dimension) {
			case DIMENSION::D3: {
				m.resize(16,0);
				m[0 *4 + 0] = x;
				m[1 *4 + 1] = y;
				m[2 *4 + 2] = z;
				m[3 *4 + 3] = 1;

			} break;

			case DIMENSION::D2: {
				m.resize(9,0);
				m[0 *3 + 0] = x;
				m[1 *3 + 1] = y;
				m[2 *3 + 2] = 1;

			} break;
			default:
				ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: unsupported operation.";
		}
}
void espreso::CoordinateSystemConfiguration::createTranslationMatrixToCenter(std::vector<double> &m) const {

	switch (dimension) {
		case DIMENSION::D3: {
			double x,y,z;

			center.x.evaluator->evaluate(1, NULL, NULL, 0, &x);
			center.y.evaluator->evaluate(1, NULL, NULL, 0, &y);
			center.z.evaluator->evaluate(1, NULL, NULL, 0, &z);

			createTranslationMatrix(m,-x,-y,-z);

		} break;

		case DIMENSION::D2: {
			double x,y;

			center.x.evaluator->evaluate(1, NULL, NULL, 0, &x);
			center.y.evaluator->evaluate(1, NULL, NULL, 0, &y);

			createTranslationMatrix(m,x,y,0);
		} break;
		default:
			ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: unsupported operation.";
	}
}


void espreso::CoordinateSystemConfiguration::createTranslationMatrixToZero(std::vector<double> &m) const{

	switch (dimension) {
		case DIMENSION::D3: {
			double x,y,z;

			center.x.evaluator->evaluate(1, NULL, NULL, 0, &x);
			center.y.evaluator->evaluate(1, NULL, NULL, 0, &y);
			center.z.evaluator->evaluate(1, NULL, NULL, 0, &z);

			createTranslationMatrix(m,x,y,z);
		} break;

		case DIMENSION::D2: {
			double x,y;

			center.x.evaluator->evaluate(1, NULL, NULL, 0, &x);
			center.y.evaluator->evaluate(1, NULL, NULL, 0, &y);

			createTranslationMatrix(m,x,y,0);
		} break;
		default:
			ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: unsupported operation.";
	}
}


void matMatProduct(
			eslocal aRows, eslocal aCols, double* aVals,
			eslocal bCols, double* bVals,
			double* cVals){
	double alpha = 1.0;
	double beta = 0.0;

	cblas_dgemm(CblasRowMajor, CblasNoTrans,
			CblasNoTrans, aRows, bCols, aCols, alpha,
			aVals, aRows, bVals, aCols, beta, cVals, aRows);
}

void espreso::CoordinateSystemConfiguration::multiplyTransformationMatrices(std::vector<double> &left, std::vector<double> &result) const {
	std::vector<double> right;
	right = result;
	result.clear();
	switch (dimension) {
		case DIMENSION::D3: {
			result.resize(16);
			//espreso::MATH::matMatProduct(4,4,left.data(),4,right.data(),tmp.data());
			matMatProduct(4,4,left.data(),4,right.data(),result.data());
		} break;
		case DIMENSION::D2: {
			result.resize(9);
			matMatProduct(3,3,left.data(),3,right.data(),result.data());
		} break;
		default:
			ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: unsupported operation.";
	}
}

espreso::Point espreso::CoordinateSystemConfiguration::applyTransformation(std::vector<double> &m, const Point &p) const{
	Point result;

	switch (dimension) {
		case DIMENSION::D3: {
			result.x = m[0*4 + 0] * p.x + m[0*4 + 1] * p.y + m[0*4 + 2] * p.z + m[0*4 + 3];
			result.y = m[1*4 + 0] * p.x + m[1*4 + 1] * p.y + m[1*4 + 2] * p.z + m[1*4 + 3];
			result.z = m[2*4 + 0] * p.x + m[2*4 + 1] * p.y + m[2*4 + 2] * p.z + m[2*4 + 3];

			result/=(m[3*4 + 0] * p.x + m[3*4 + 1] * p.y + m[3*4 + 2] * p.z + m[3*4 + 3]);
		} break;
		case DIMENSION::D2: {
			result.x = m[0*3 + 0] * p.x + m[0*3 + 1] * p.y + m[0*3 + 2];
			result.y = m[1*3 + 0] * p.x + m[1*3 + 1] * p.y + m[1*3 + 2];

			result/=(m[2*3 + 0] * p.x + m[2*3 + 1] * p.y + m[2*3 + 2]);
		} break;
		default:
			ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: unsupported operation.";
	}
	return result;
}


void  espreso::CoordinateSystemConfiguration::createRotationMatrix(std::vector<double> &m) const {
	auto d2r = [] (double degree) -> double {
			return M_PI * degree / 180;
	};
	Point rPoint;
	switch (dimension) {
		case DIMENSION::D3: {
			rotation.x.evaluator->evaluate(1, NULL, NULL, 0, &(rPoint.x));
			rotation.y.evaluator->evaluate(1, NULL, NULL, 0, &(rPoint.y));
			rotation.z.evaluator->evaluate(1, NULL, NULL, 0, &(rPoint.z));

		} break;
		case DIMENSION::Z:
		case DIMENSION::D2: {
			rotation.z.evaluator->evaluate(1, NULL, NULL, 0, &(rPoint.z));
		} break;
		default:
			ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: unsupported operation.";
	}

	switch (dimension) {
		case DIMENSION::D3: {
			m.clear();
			m.resize(16,0);
			Point sin, cos;
			switch (type) {
				case CoordinateSystemConfiguration::TYPE::CARTESIAN:{
					cos.x = std::cos(d2r(rPoint.x));
					cos.y = std::cos(d2r(rPoint.y));
					cos.z = std::cos(d2r(rPoint.z));

					sin.x = std::sin(d2r(rPoint.x));
					sin.y = std::sin(d2r(rPoint.y));
					sin.z = std::sin(d2r(rPoint.z));
				}break;
				default:
					ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: unsupported operation.";

			}

			m[0*4 + 0] = cos.y * cos.z;                      	m[0*4 + 1] = cos.y * sin.z;                         m[0*4 + 2] = -sin.y;
			m[1*4 + 0] = cos.z * sin.x * sin.y - cos.x * sin.z; m[1*4 + 1] = cos.x * cos.z + sin.x * sin.y * sin.z; m[1*4 + 2] = cos.y * sin.x;
			m[2*4 + 0] = sin.x * sin.z + cos.x * cos.z * sin.y; m[2*4 + 1] = cos.x * sin.y * sin.z - cos.z * sin.x; m[2*4 + 2] = cos.x * cos.y;
			m[3*4 + 3] = 1;

		} break;

		case DIMENSION::Z:
		case DIMENSION::D2: {
			double cos, sin;
			m.clear();
			m.resize(9,0);

			switch (type) {

				case CoordinateSystemConfiguration::TYPE::CARTESIAN: {
					cos = std::cos(d2r(rPoint.z));
					sin = std::sin(d2r(rPoint.z));
				}break;

				default:
					ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: unsupported operation.";
			}

			m[0*3 + 0] =  cos; m[0*3 + 1] = sin;
			m[1*3 + 0] = -sin; m[1*3 + 1] = cos;
			m[2*3 + 2] = 1;
		}
		default:
			ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: unsupported operation.";
	}
}




