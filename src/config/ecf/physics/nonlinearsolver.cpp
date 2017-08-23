
#include "nonlinearsolver.h"
#include "../../configuration.hpp"

espreso::NonLinearSolverConfiguration::NonLinearSolverConfiguration(const std::string &solution_name, const std::string &residual_name)
{
	method = METHOD::NEWTON_RHAPSON;
	REGISTER(method, ECFMetaData()
			.setdescription({ "Non-linear method." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("NEWTON_RHAPSON").setdescription("Newton-Rhapson."))
			.addoption(ECFOption().setname("MODIFIED_NEWTON_RHAPSON").setdescription("Newton-Rhapson without re-assembling of stiffness matrices.")));

	stepping = STEPPINGG::FALSE;
	REGISTER(stepping, ECFMetaData()
			.setdescription({ "Use substeps." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("TRUE").setdescription("Turn on."))
			.addoption(ECFOption().setname("FALSE").setdescription("Turn off."))
			.addoption(ECFOption().setname("AUTO").setdescription("Automatic.")));

	max_iterations = 15;
	REGISTER(max_iterations, ECFMetaData()
			.setdescription({ "Maximal number on non-linear iterations." })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	substeps = 1;
	REGISTER(substeps, ECFMetaData()
			.setdescription({ "Number of substeps when stepping is true." })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

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

	check_solution = true;
	check_residual = false;
	registerParameter("check_" + solution_name, check_solution, ECFMetaData()
			.setdescription({ "Stopping criteria based on solution." })
			.setdatatype({ ECFDataType::FLOAT }));
	registerParameter("check_" + residual_name, check_residual, ECFMetaData()
			.setdescription({ "Stopping criteria based on residual." })
			.setdatatype({ ECFDataType::FLOAT }));

	requested_solution = requested_residual = 1e-3;
	registerParameter("requested_" + solution_name, requested_solution, ECFMetaData()
			.setdescription({ "Requested solution precision." })
			.setdatatype({ ECFDataType::FLOAT }));
	registerParameter("requested_" + residual_name, requested_residual, ECFMetaData()
			.setdescription({ "Requested residual precision." })
			.setdatatype({ ECFDataType::FLOAT }));

	r_tol = 0.1;
	c_fact = 0.8;
	REGISTER(r_tol, ECFMetaData()
			.setdescription({ "Parameter for decrease outer precision." })
			.setdatatype({ ECFDataType::FLOAT }));
	REGISTER(c_fact, ECFMetaData()
			.setdescription({ "Parameter for decrease inner precision." })
			.setdatatype({ ECFDataType::FLOAT }));
}



