
#include "transientsecondorderimplicit.h"
#include "config/configuration.hpp"

using namespace espreso;

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

	numerical_dumping = 0;
	REGISTER(numerical_dumping, ECFMetaData()
			.setdescription({ "Numerical Dumping" })
			.setdatatype({ ECFDataType::FLOAT }));

	ecfdescription->addSpace();

	dumping = DUMPING::NONE;
	REGISTER(method, ECFMetaData()
			.setdescription({ "Dumping type" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("NONE").setdescription("Without dumping."))
			.addoption(ECFOption().setname("DIRECT").setdescription("Direct."))
			.addoption(ECFOption().setname("DUMPING_RATIO").setdescription("By ratio.")));

	stiffness_dumping = 0;
	REGISTER(stiffness_dumping, ECFMetaData()
			.setdescription({ "Stiffnes Dumping" })
			.setdatatype({ ECFDataType::FLOAT }));

	mass_dumping = 0;
	REGISTER(mass_dumping, ECFMetaData()
			.setdescription({ "Mass Dumping" })
			.setdatatype({ ECFDataType::FLOAT }));

	dumping_ratio = 0;
	REGISTER(dumping_ratio, ECFMetaData()
			.setdescription({ "Dumping Ratio" })
			.setdatatype({ ECFDataType::FLOAT }));

	frequency = 0;
	REGISTER(frequency, ECFMetaData()
			.setdescription({ "Dumping frequency" })
			.setdatatype({ ECFDataType::FLOAT }));

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
