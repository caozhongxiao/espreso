
#include "coordinatesystem.h"

#include "../../configuration.hpp"

espreso::CoordinateSystemConfiguration::CoordinateSystemConfiguration()
{
	type = TYPE::CARTESIAN;
	dimension = DIMENSION::D3;
	rotation_x.value = rotation_y.value = rotation_z.value = center_x.value = center_y.value = center_z.value = "0";

	REGISTER(type, ECFMetaData()
			.setdescription({ "Type" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("CARTESIAN").setdescription("Cartesian system."))
			.addoption(ECFOption().setname("CYLINDRICAL").setdescription("Cylindrical system."))
			.addoption(ECFOption().setname("SPHERICAL").setdescription("Spherical system.").allowonly([&] () { return dimension == DIMENSION::D3; })));

	REGISTER(rotation_x, ECFMetaData()
			.setdescription({ "A x-rotation of the material." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setcoordinatevariables()
			.allowonly([&] () { return type == TYPE::CARTESIAN && dimension == DIMENSION::D3; }));
	REGISTER(rotation_y, ECFMetaData()
			.setdescription({ "A y-rotation of the material." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setcoordinatevariables()
			.allowonly([&] () { return type == TYPE::CARTESIAN && dimension == DIMENSION::D3; }));
	REGISTER(rotation_z, ECFMetaData()
			.setdescription({ "A z-rotation of the material." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setcoordinatevariables()
			.allowonly([&] () { return type == TYPE::CARTESIAN; }));

	REGISTER(center_x, ECFMetaData()
			.setdescription({ "A x-center of the material." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setcoordinatevariables()
			.allowonly([&] () { return type != TYPE::CARTESIAN; }));
	REGISTER(center_y, ECFMetaData()
			.setdescription({ "A y-center of the material." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setcoordinatevariables()
			.allowonly([&] () { return type != TYPE::CARTESIAN; }));
	REGISTER(center_z, ECFMetaData()
			.setdescription({ "A z-center of the material." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setcoordinatevariables()
			.allowonly([&] () { return type == TYPE::SPHERICAL; }));
}



