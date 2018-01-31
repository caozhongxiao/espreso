
#include "coordinatesystem.h"

#include "../../configuration.hpp"

espreso::CoordinateSystemConfiguration::CoordinateSystemConfiguration()
: dimension(DIMENSION::D3), rotation(dimension, true), center(dimension, true)
{
	type = TYPE::CARTESIAN;

	REGISTER(type, ECFMetaData()
			.setdescription({ "Type" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("CARTESIAN").setdescription("Cartesian system."))
			.addoption(ECFOption().setname("CYLINDRICAL").setdescription("Cylindrical system."))
			.addoption(ECFOption().setname("SPHERICAL").setdescription("Spherical system.").allowonly([&] () { return dimension == DIMENSION::D3; })));

	REGISTER(rotation, ECFMetaData()
			.setdescription({ "A x-rotation of the material." })
			.setcoordinatevariables()
			.allowonly([&] () { return type == TYPE::CARTESIAN; }));

	REGISTER(center, ECFMetaData()
			.setdescription({ "A x-center of the material." })
			.setcoordinatevariables()
			.allowonly([&] () { return type != TYPE::CARTESIAN; }));
}



