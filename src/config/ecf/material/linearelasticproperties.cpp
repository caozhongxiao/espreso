
#include "linearelasticproperties.h"
#include "../../configuration.hpp"

espreso::LinearElasticPropertiesConfiguration::LinearElasticPropertiesConfiguration()
: model(MODEL::ISOTROPIC),
  is3D(true),
  poisson_ratio(3),
  young_modulus(3),
  thermal_expansion(3),
  shear_modulus(3),
  anisotropic(6)
{
	REGISTER(model, ECFMetaData()
			.setdescription({ "Tensors model." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("ISOTROPIC").setdescription("Isotropic model."))
			.addoption(ECFOption().setname("ORTHOTROPIC").setdescription("Orthotropic model."))
			.addoption(ECFOption().setname("ANISOTROPIC").setdescription("Anisotropic model.").allowonly([&] () { return is3D; })));

	addSeparator();

	registerParameter("MIXY", poisson_ratio.value(0, 0), ECFMetaData()
			.setdescription({ "Poisson ratio XY." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(poisson_ratio)
			.allowonly([&] () { return model != MODEL::ANISOTROPIC; }));
	registerParameter("MIXZ", poisson_ratio.value(1, 1), ECFMetaData()
			.setdescription({ "Poisson ratio XZ." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(poisson_ratio)
			.allowonly([&] () { return model == MODEL::ORTHOTROPIC && is3D; }));
	registerParameter("MIYZ", poisson_ratio.value(2, 2), ECFMetaData()
			.setdescription({ "Poisson ratio YZ." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(poisson_ratio)
			.allowonly([&] () { return model == MODEL::ORTHOTROPIC && is3D; }));

	addSeparator();

	registerParameter("EX", young_modulus.value(0, 0), ECFMetaData()
			.setdescription({ "Young modulus X." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(young_modulus)
			.allowonly([&] () { return model != MODEL::ANISOTROPIC; }));
	registerParameter("EY", young_modulus.value(1, 1), ECFMetaData()
			.setdescription({ "Young modulus Y." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(young_modulus)
			.allowonly([&] () { return model == MODEL::ORTHOTROPIC; }));
	registerParameter("EZ", young_modulus.value(2, 2), ECFMetaData()
			.setdescription({ "Young modulus Z." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(young_modulus)
			.allowonly([&] () { return model == MODEL::ORTHOTROPIC && is3D; }));

	addSeparator();

	registerParameter("TEX", thermal_expansion.value(0, 0), ECFMetaData()
			.setdescription({ "Thermal expansion X." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(thermal_expansion));

	registerParameter("TEY", thermal_expansion.value(1, 1), ECFMetaData()
			.setdescription({ "Thermal expansion Y." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(thermal_expansion)
			.allowonly([&] () { return model != MODEL::ISOTROPIC; }));
	registerParameter("TEZ", thermal_expansion.value(2, 2), ECFMetaData()
			.setdescription({ "Thermal expansion Z." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(thermal_expansion)
			.allowonly([&] () { return model != MODEL::ISOTROPIC && is3D; }));

	addSeparator();

	registerParameter("GXY", shear_modulus.value(0, 0), ECFMetaData()
			.setdescription({ "Shear modulus XY." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(shear_modulus)
			.allowonly([&] () { return model == MODEL::ORTHOTROPIC && is3D; }));
	registerParameter("GXZ", shear_modulus.value(1, 1), ECFMetaData()
			.setdescription({ "Shear modulus XZ." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(shear_modulus)
			.allowonly([&] () { return model == MODEL::ORTHOTROPIC && is3D; }));
	registerParameter("GYZ", shear_modulus.value(2, 2), ECFMetaData()
			.setdescription({ "Shear modulus YZ." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(shear_modulus)
			.allowonly([&] () { return model == MODEL::ORTHOTROPIC && is3D; }));

	addSeparator();

	// anisotropic is allowed only in 3D
	registerParameter("D11", anisotropic.value(0, 0), ECFMetaData()
			.setdescription({ "Anisotropic parameter (1, 1)." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(anisotropic)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC; }));
	registerParameter("D12", anisotropic.value(0, 1), ECFMetaData()
			.setdescription({ "Anisotropic parameter (1, 2)." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(anisotropic)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC; }));
	registerParameter("D13", anisotropic.value(0, 2), ECFMetaData()
			.setdescription({ "Anisotropic parameter (1, 3)." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(anisotropic)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC; }));
	registerParameter("D14", anisotropic.value(0, 3), ECFMetaData()
			.setdescription({ "Anisotropic parameter (1, 4)." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(anisotropic)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC; }));
	registerParameter("D15", anisotropic.value(0, 4), ECFMetaData()
			.setdescription({ "Anisotropic parameter (1, 5)." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(anisotropic)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC; }));
	registerParameter("D16", anisotropic.value(0, 5), ECFMetaData()
			.setdescription({ "Anisotropic parameter (1, 6)." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(anisotropic)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC; }));
	registerParameter("D22", anisotropic.value(1, 1), ECFMetaData()
			.setdescription({ "Anisotropic parameter (2, 2)." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(anisotropic)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC; }));
	registerParameter("D23", anisotropic.value(1, 2), ECFMetaData()
			.setdescription({ "Anisotropic parameter (2, 3)." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(anisotropic)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC; }));
	registerParameter("D24", anisotropic.value(1, 3), ECFMetaData()
			.setdescription({ "Anisotropic parameter (2, 4)." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(anisotropic)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC; }));
	registerParameter("D25", anisotropic.value(1, 4), ECFMetaData()
			.setdescription({ "Anisotropic parameter (2, 5)." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(anisotropic)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC; }));
	registerParameter("D26", anisotropic.value(1, 5), ECFMetaData()
			.setdescription({ "Anisotropic parameter (2, 6)." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(anisotropic)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC; }));
	registerParameter("D33", anisotropic.value(2, 2), ECFMetaData()
			.setdescription({ "Anisotropic parameter (3, 3)." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(anisotropic)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC; }));
	registerParameter("D34", anisotropic.value(2, 3), ECFMetaData()
			.setdescription({ "Anisotropic parameter (3, 4)." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(anisotropic)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC; }));
	registerParameter("D35", anisotropic.value(2, 4), ECFMetaData()
			.setdescription({ "Anisotropic parameter (3, 5)." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(anisotropic)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC; }));
	registerParameter("D36", anisotropic.value(2, 5), ECFMetaData()
			.setdescription({ "Anisotropic parameter (3, 6)." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(anisotropic)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC; }));
	registerParameter("D44", anisotropic.value(3, 3), ECFMetaData()
			.setdescription({ "Anisotropic parameter (4, 4)." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(anisotropic)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC; }));
	registerParameter("D45", anisotropic.value(3, 4), ECFMetaData()
			.setdescription({ "Anisotropic parameter (4, 5)." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(anisotropic)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC; }));
	registerParameter("D46", anisotropic.value(3, 5), ECFMetaData()
			.setdescription({ "Anisotropic parameter (4, 6)." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(anisotropic)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC; }));
	registerParameter("D55", anisotropic.value(4, 4), ECFMetaData()
			.setdescription({ "Anisotropic parameter (5, 5)." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(anisotropic)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC; }));
	registerParameter("D56", anisotropic.value(4, 5), ECFMetaData()
			.setdescription({ "Anisotropic parameter (5, 6)." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(anisotropic)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC; }));
	registerParameter("D66", anisotropic.value(5, 5), ECFMetaData()
			.setdescription({ "Anisotropic parameter (6, 6)." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setmaterialvariables()
			.settensor(anisotropic)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC; }));
}

