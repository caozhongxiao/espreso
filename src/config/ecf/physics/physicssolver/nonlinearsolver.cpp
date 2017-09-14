
#include "nonlinearsolver.h"
#include "../../../configuration.hpp"

espreso::NonLinearSolverConfiguration::NonLinearSolverConfiguration(const std::string &firstResidualName, const std::string &secondResidualName)
{
	method = METHOD::NEWTON_RAPHSON;
	REGISTER(method, ECFMetaData()
			.setdescription({ "Non-linear method." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("NEWTON_RAPHSON").setdescription("Newton-Raphson."))
			.addoption(ECFOption().setname("MODIFIED_NEWTON_RAPHSON").setdescription("Newton-Raphson without re-assembling of stiffness matrices.")));

	addSpace();

	check_first_residual = true;
	check_second_residual = false;
	registerParameter("check_" + firstResidualName, check_first_residual, ECFMetaData()
			.setdescription({ "Stopping criteria based on solution." })
			.setdatatype({ ECFDataType::FLOAT }));
	registerParameter("check_" + secondResidualName, check_second_residual, ECFMetaData()
			.setdescription({ "Stopping criteria based on residual." })
			.setdatatype({ ECFDataType::FLOAT }));

	requested_first_residual = requested_second_residual = 1e-3;
	registerParameter("requested_" + firstResidualName + "_residual", requested_first_residual, ECFMetaData()
			.setdescription({ "Requested solution precision." })
			.setdatatype({ ECFDataType::FLOAT }));
	registerParameter("requested_" + secondResidualName + "_residual", requested_second_residual, ECFMetaData()
			.setdescription({ "Requested residual precision." })
			.setdatatype({ ECFDataType::FLOAT }));

	addSeparator();

	stepping = STEPPINGG::FALSE;
	REGISTER(stepping, ECFMetaData()
			.setdescription({ "Use substeps." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("TRUE").setdescription("Turn on."))
			.addoption(ECFOption().setname("FALSE").setdescription("Turn off."))
			.addoption(ECFOption().setname("AUTO").setdescription("Automatic.")));

	substeps = 1;
	REGISTER(substeps, ECFMetaData()
			.setdescription({ "Number of substeps when stepping is true." })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	max_iterations = 15;
	REGISTER(max_iterations, ECFMetaData()
			.setdescription({ "Maximal number on non-linear iterations." })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	addSpace();

	line_search = tangent_matrix_correction = adaptive_precision = false;
	REGISTER(line_search, ECFMetaData()
			.setdescription({ "Use line search." })
			.setdatatype({ ECFDataType::BOOL }));
	REGISTER(tangent_matrix_correction, ECFMetaData()
			.setdescription({ "Correction of stiffness matrix." })
			.setdatatype({ ECFDataType::BOOL }));
	REGISTER(adaptive_precision, ECFMetaData()
			.setdescription({ "Adaptive precision of a linear solver." })
			.setdatatype({ ECFDataType::BOOL }));

	addSpace();

	r_tol = 0.1;
	c_fact = 0.8;
	REGISTER(r_tol, ECFMetaData()
			.setdescription({ "Parameter for decrease outer precision." })
			.setdatatype({ ECFDataType::FLOAT }));
	REGISTER(c_fact, ECFMetaData()
			.setdescription({ "Parameter for decrease inner precision." })
			.setdatatype({ ECFDataType::FLOAT }));
}