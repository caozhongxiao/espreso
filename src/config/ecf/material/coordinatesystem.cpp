
#include "coordinatesystem.h"

#include "../../configuration.hpp"

espreso::CoordinateSystemConfiguration::CoordinateSystemConfiguration()
{
	type = TYPE::CARTESIAN;

	REGISTER(rotation_x, ECFMetaData()
			.setdescription({ "A x-rotation of the material." })
			.setdatatype({ ECFDataType::EXPRESSION }));
	REGISTER(rotation_y, ECFMetaData()
			.setdescription({ "A y-rotation of the material." })
			.setdatatype({ ECFDataType::EXPRESSION }));
	REGISTER(rotation_z, ECFMetaData()
			.setdescription({ "A z-rotation of the material." })
			.setdatatype({ ECFDataType::EXPRESSION }));

	REGISTER(center_x, ECFMetaData()
			.setdescription({ "A x-center of the material." })
			.setdatatype({ ECFDataType::EXPRESSION }));
	REGISTER(center_y, ECFMetaData()
			.setdescription({ "A y-center of the material." })
			.setdatatype({ ECFDataType::EXPRESSION }));
	REGISTER(center_z, ECFMetaData()
			.setdescription({ "A z-center of the material." })
			.setdatatype({ ECFDataType::EXPRESSION }));
}



