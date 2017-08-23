
#include "thermalproperties.h"
#include "../../configuration.hpp"

espreso::ThermalPropertiesConfiguration::ThermalPropertiesConfiguration()
: model(MODEL::ISOTROPIC),
  is3D(true),
  thermal_conductivity(3)
{
	REGISTER(model, ECFMetaData()
			.setdescription({ "Tensors model." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("ISOTROPIC").setdescription("Isotropic model."))
			.addoption(ECFOption().setname("DIAGONAL").setdescription("Orthotropic model."))
			.addoption(ECFOption().setname("SYMMETRIC").setdescription("Anisotropic model."))
			.addoption(ECFOption().setname("ANISOTROPIC").setdescription("Anisotropic model.").allowonly([&] () { return is3D; })));

	addSeparator();

	registerParameter("KXX", thermal_conductivity.value(0, 0), ECFMetaData()
			.setdescription({ "Thermal conductivity XX." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(thermal_conductivity));
	registerParameter("KYY", thermal_conductivity.value(1, 1), ECFMetaData()
			.setdescription({ "Thermal conductivity YY." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(thermal_conductivity)
			.allowonly([&] () { return model != MODEL::ISOTROPIC; }));
	registerParameter("KZZ", thermal_conductivity.value(2, 2), ECFMetaData()
			.setdescription({ "Thermal conductivity ZZ." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(thermal_conductivity)
			.allowonly([&] () { return model != MODEL::ISOTROPIC && is3D; }));

	registerParameter("KXY", thermal_conductivity.value(0, 1), ECFMetaData()
			.setdescription({ "Thermal conductivity XY." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(thermal_conductivity)
			.allowonly([&] () { return model == MODEL::SYMMETRIC || model == MODEL::ANISOTROPIC; }));
	registerParameter("KXZ", thermal_conductivity.value(0, 2), ECFMetaData()
			.setdescription({ "Thermal conductivity XZ." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(thermal_conductivity)
			.allowonly([&] () { return (model == MODEL::SYMMETRIC || model == MODEL::ANISOTROPIC) && is3D; }));
	registerParameter("KYZ", thermal_conductivity.value(1, 2), ECFMetaData()
			.setdescription({ "Thermal conductivity YZ." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(thermal_conductivity)
			.allowonly([&] () { return (model == MODEL::SYMMETRIC || model == MODEL::ANISOTROPIC)  && is3D; }));

	registerParameter("KYX", thermal_conductivity.value(1, 0), ECFMetaData()
			.setdescription({ "Thermal conductivity YX." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(thermal_conductivity)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC; }));
	registerParameter("KZX", thermal_conductivity.value(2, 0), ECFMetaData()
			.setdescription({ "Thermal conductivity ZX." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(thermal_conductivity)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC && is3D; }));
	registerParameter("KZY", thermal_conductivity.value(2, 1), ECFMetaData()
			.setdescription({ "Thermal conductivity ZY." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(thermal_conductivity)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC && is3D; }));
}
