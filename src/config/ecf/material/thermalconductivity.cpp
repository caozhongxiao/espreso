
#include "config/configuration.hpp"
#include "config/conditions.h"
#include "thermalconductivity.h"

espreso::ThermalConductivityConfiguration::ThermalConductivityConfiguration(DIMENSION *D)
: model(MODEL::ISOTROPIC),
  dimension(D),
  values(3)
{
	REGISTER(model, ECFMetaData()
			.setdescription({ "Model" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("ISOTROPIC").setdescription("Isotropic model."))
			.addoption(ECFOption().setname("DIAGONAL").setdescription("Orthotropic model."))
			.addoption(ECFOption().setname("SYMMETRIC").setdescription("Symmetric model."))
			.addoption(ECFOption().setname("ANISOTROPIC").setdescription("Anisotropic model.").allowonly([&] () { return *dimension == DIMENSION::D3; })));

	ecfdescription->addSeparator();

	ecfdescription->registerParameter("KXX", values.get(0, 0), ECFMetaData()
			.setdescription({ "KXX" })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(values));
	ecfdescription->registerParameter("KYY", values.get(1, 1), ECFMetaData()
			.setdescription({ "KYY" })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(values)
			.allowonly([&] () { return model != MODEL::ISOTROPIC; })
			.addconstraint(!ECFCondition(model, MODEL::ISOTROPIC)));
	ecfdescription->registerParameter("KZZ", values.get(2, 2), ECFMetaData()
			.setdescription({ "KZZ" })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(values)
			.allowonly([&] () { return model != MODEL::ISOTROPIC && *dimension == DIMENSION::D3; }));

	ecfdescription->registerParameter("KXY", values.get(0, 1), ECFMetaData()
			.setdescription({ "KXY" })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(values)
			.allowonly([&] () { return model == MODEL::SYMMETRIC || model == MODEL::ANISOTROPIC; }));
	ecfdescription->registerParameter("KXZ", values.get(0, 2), ECFMetaData()
			.setdescription({ "KXZ" })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(values)
			.allowonly([&] () { return (model == MODEL::SYMMETRIC || model == MODEL::ANISOTROPIC) && *dimension == DIMENSION::D3; }));
	ecfdescription->registerParameter("KYZ", values.get(1, 2), ECFMetaData()
			.setdescription({ "KYZ" })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(values)
			.allowonly([&] () { return (model == MODEL::SYMMETRIC || model == MODEL::ANISOTROPIC)  && *dimension == DIMENSION::D3; }));

	ecfdescription->registerParameter("KYX", values.get(1, 0), ECFMetaData()
			.setdescription({ "KYX" })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(values)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC; }));
	ecfdescription->registerParameter("KZX", values.get(2, 0), ECFMetaData()
			.setdescription({ "KZX" })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(values)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC && *dimension == DIMENSION::D3; }));
	ecfdescription->registerParameter("KZY", values.get(2, 1), ECFMetaData()
			.setdescription({ "KZY" })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(values)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC && *dimension == DIMENSION::D3; }));
}
