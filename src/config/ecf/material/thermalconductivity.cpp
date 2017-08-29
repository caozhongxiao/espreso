
#include "../../configuration.hpp"
#include "thermalconductivity.h"

espreso::ThermalConductivityConfiguration::ThermalConductivityConfiguration()
: model(MODEL::ISOTROPIC),
  is3D(true),
  values(3)
{
	REGISTER(model, ECFMetaData()
			.setdescription({ "Tensors model." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("ISOTROPIC").setdescription("Isotropic model."))
			.addoption(ECFOption().setname("DIAGONAL").setdescription("Orthotropic model."))
            .addoption(ECFOption().setname("SYMMETRIC").setdescription("Symmetric model."))
			.addoption(ECFOption().setname("ANISOTROPIC").setdescription("Anisotropic model.").allowonly([&] () { return is3D; })));

	addSeparator();

	registerParameter("KXX", values.get(0, 0), ECFMetaData()
			.setdescription({ "Thermal conductivity XX." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(values));
	registerParameter("KYY", values.get(1, 1), ECFMetaData()
			.setdescription({ "Thermal conductivity YY." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(values)
			.allowonly([&] () { return model != MODEL::ISOTROPIC; }));
	registerParameter("KZZ", values.get(2, 2), ECFMetaData()
			.setdescription({ "Thermal conductivity ZZ." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(values)
			.allowonly([&] () { return model != MODEL::ISOTROPIC && is3D; }));

	registerParameter("KXY", values.get(0, 1), ECFMetaData()
			.setdescription({ "Thermal conductivity XY." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(values)
			.allowonly([&] () { return model == MODEL::SYMMETRIC || model == MODEL::ANISOTROPIC; }));
	registerParameter("KXZ", values.get(0, 2), ECFMetaData()
			.setdescription({ "Thermal conductivity XZ." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(values)
			.allowonly([&] () { return (model == MODEL::SYMMETRIC || model == MODEL::ANISOTROPIC) && is3D; }));
	registerParameter("KYZ", values.get(1, 2), ECFMetaData()
			.setdescription({ "Thermal conductivity YZ." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(values)
			.allowonly([&] () { return (model == MODEL::SYMMETRIC || model == MODEL::ANISOTROPIC)  && is3D; }));

	registerParameter("KYX", values.get(1, 0), ECFMetaData()
			.setdescription({ "Thermal conductivity YX." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(values)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC; }));
	registerParameter("KZX", values.get(2, 0), ECFMetaData()
			.setdescription({ "Thermal conductivity ZX." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(values)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC && is3D; }));
	registerParameter("KZY", values.get(2, 1), ECFMetaData()
			.setdescription({ "Thermal conductivity ZY." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(values)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC && is3D; }));
}
