
#include "gridtower.h"

#include "../../configuration.hpp"

espreso::GridTowerGeneratorConfiguration::GridTowerGeneratorConfiguration()
{
	direction = DIRECTION::X;

	grids[0] = GridGeneratorConfiguration();

	REGISTER(direction, ECFMetaData()
			.setdescription({ "Direction of generated tower." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("X").setdescription("Grids are placed in x-direction."))
			.addoption(ECFOption().setname("Y").setdescription("Grids are placed in y-direction."))
			.addoption(ECFOption().setname("Z").setdescription("Grids are placed in z-direction.")));

	REGISTER(grids, ECFMetaData()
			.setdescription({ "An index of grid in tower. Indices has to be continuous starting from 0.", "Description of grid." })
			.setdatatype({ ECFDataType::NONNEGATIVE_INTEGER })
			.setpattern({ "0" }));
}


