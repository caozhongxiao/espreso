
#include "transientsecondorderimplicit.h"
#include "config/configuration.hpp"

using namespace espreso;

TransientSecondOrderImplicitDirectDampingConfiguration::TransientSecondOrderImplicitDirectDampingConfiguration()
: mass({ "TIME" }), stiffness({ "TIME" })
{
	stiffness.value = "0";
	REGISTER(stiffness, ECFMetaData()
			.setdescription({ "Stiffnes Damping" })
			.setdatatype({ ECFDataType::FLOAT }));

	mass.value = "0";
	REGISTER(mass, ECFMetaData()
			.setdescription({ "Mass Damping" })
			.setdatatype({ ECFDataType::FLOAT }));
}

TransientSecondOrderImplicitRatioDampingConfiguration::TransientSecondOrderImplicitRatioDampingConfiguration()
: ratio({ "TIME" }), frequency({ "TIME" })
{
	ratio.value = "0";
	REGISTER(ratio, ECFMetaData()
			.setdescription({ "Damping Ratio" })
			.setdatatype({ ECFDataType::FLOAT }));

	frequency.value = "0";
	REGISTER(frequency, ECFMetaData()
			.setdescription({ "Damping frequency" })
			.setdatatype({ ECFDataType::FLOAT }));
}

TransientSecondOrderImplicitConfiguration::TransientSecondOrderImplicitConfiguration()
{
	method = METHOD::NEWMARK;
	REGISTER(method, ECFMetaData()
			.setdescription({ "Method" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("NEWMARK").setdescription("NEWMARK method.")));

	alpha = 0.5;
	REGISTER(alpha, ECFMetaData()
			.setdescription({ "Alpha" })
			.setdatatype({ ECFDataType::FLOAT }));

	delta = 0.25;
	REGISTER(delta, ECFMetaData()
			.setdescription({ "Delta" })
			.setdatatype({ ECFDataType::FLOAT }));

	numerical_damping = 0;
	REGISTER(numerical_damping, ECFMetaData()
			.setdescription({ "Numerical Damping" })
			.setdatatype({ ECFDataType::FLOAT }));

	ecfdescription->addSpace();

	damping = DAMPING::NONE;
	REGISTER(damping, ECFMetaData()
			.setdescription({ "Damping type" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("NONE").setdescription("Without damping."))
			.addoption(ECFOption().setname("DIRECT").setdescription("Direct."))
			.addoption(ECFOption().setname("DAMPING_RATIO").setdescription("By ratio.")));

	REGISTER(direct_damping, ECFMetaData()
			.setdescription({ "Direct damping stepping" }));

	REGISTER(ratio_damping, ECFMetaData()
			.setdescription({ "Ratio damping stepping" }));

	mass_matrix_type = MASS_MATRIX_TYPE::CONSISTENT;
	REGISTER(method, ECFMetaData()
			.setdescription({ "Method" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("CONSISTENT").setdescription("Consistent method."))
			.addoption(ECFOption().setname("DIAGONAL").setdescription("Diagonal method."))
			.addoption(ECFOption().setname("HRZDIAGONAL").setdescription("HRZ Diagonal method.")));

	ecfdescription->addSpace();

	time_step = 0.1;
	REGISTER(time_step, ECFMetaData()
			.setdescription({ "Time step" })
			.setdatatype({ ECFDataType::FLOAT }));

	REGISTER(auto_time_stepping, ECFMetaData()
			.setdescription({ "Auto time stepping" }));
}
