
#include "hyperelasticproperties.h"
#include "config/configuration.hpp"

espreso::HyperElasticPropertiesConfiguration::HyperElasticPropertiesConfiguration()
: model(MODEL::NEO_HOOKEN),
  dimension(DIMENSION::D3),
  mu(ECFMetaData::getmaterialvariables()),
  d(ECFMetaData::getmaterialvariables()),
  C10(ECFMetaData::getmaterialvariables()),
  C01(ECFMetaData::getmaterialvariables()),
  C11(ECFMetaData::getmaterialvariables()),
  C02(ECFMetaData::getmaterialvariables()),
  C20(ECFMetaData::getmaterialvariables()),
  C30(ECFMetaData::getmaterialvariables()),
  C21(ECFMetaData::getmaterialvariables()),
  C12(ECFMetaData::getmaterialvariables()),
  C03(ECFMetaData::getmaterialvariables())
{
	REGISTER(model, ECFMetaData()
			.setdescription({ "Material model." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("NEO_HOOKEN").setdescription("Neo-Hooken"))
			.addoption(ECFOption().setname("MOONEY_RIVLIN_2PARAMS").setdescription("Mooney-Rivlin with 2 parameters"))
			.addoption(ECFOption().setname("MOONEY_RIVLIN_3PARAMS").setdescription("Mooney-Rivlin with 3 parameters"))
			.addoption(ECFOption().setname("MOONEY_RIVLIN_5PARAMS").setdescription("Mooney-Rivlin with 5 parameters"))
			.addoption(ECFOption().setname("MOONEY_RIVLIN_9PARAMS").setdescription("Mooney-Rivlin with 9 parameters")));

	ecfdescription->addSeparator();

	REGISTER(mu, ECFMetaData()
			.setdescription({ "Initial shear modulus." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.allowonly([&] () { return model == MODEL::NEO_HOOKEN; }));

	auto isMR2 = [&] () {
		return
				model == MODEL::MOONEY_RIVLIN_2PARAMS ||
				model == MODEL::MOONEY_RIVLIN_3PARAMS ||
				model == MODEL::MOONEY_RIVLIN_5PARAMS ||
				model == MODEL::MOONEY_RIVLIN_9PARAMS;
	};

	auto isMR3 = [&] () {
		return
				model == MODEL::MOONEY_RIVLIN_3PARAMS ||
				model == MODEL::MOONEY_RIVLIN_5PARAMS ||
				model == MODEL::MOONEY_RIVLIN_9PARAMS;
	};

	auto isMR5 = [&] () {
		return
				model == MODEL::MOONEY_RIVLIN_5PARAMS ||
				model == MODEL::MOONEY_RIVLIN_9PARAMS;
	};

	auto isMR9 = [&] () {
		return
				model == MODEL::MOONEY_RIVLIN_9PARAMS;
	};

	REGISTER(C01, ECFMetaData()
			.setdescription({ "Material constant." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.allowonly([&] () { return isMR2(); }));

	REGISTER(C01, ECFMetaData()
			.setdescription({ "Material constant." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.allowonly([&] () { return isMR2(); }));

	REGISTER(C11, ECFMetaData()
			.setdescription({ "Material constant." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.allowonly([&] () { return isMR3(); }));

	REGISTER(C02, ECFMetaData()
			.setdescription({ "Material constant." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.allowonly([&] () { return isMR5(); }));

	REGISTER(C20, ECFMetaData()
			.setdescription({ "Material constant." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.allowonly([&] () { return isMR5(); }));

	REGISTER(C30, ECFMetaData()
			.setdescription({ "Material constant." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.allowonly([&] () { return isMR9(); }));

	REGISTER(C21, ECFMetaData()
			.setdescription({ "Material constant." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.allowonly([&] () { return isMR9(); }));

	REGISTER(C12, ECFMetaData()
			.setdescription({ "Material constant." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.allowonly([&] () { return isMR9(); }));

	REGISTER(C03, ECFMetaData()
			.setdescription({ "Material constant." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.allowonly([&] () { return isMR9(); }));

	REGISTER(d, ECFMetaData()
			.setdescription({ "Incompressibility parameter." })
			.setdatatype({ ECFDataType::EXPRESSION }));
}

