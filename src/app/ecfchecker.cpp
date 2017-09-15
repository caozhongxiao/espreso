
#include "mpi.h"

#include "../config/ecf/ecf.h"
#include "../config/reader/reader.h"
#include "../basis/logging/logging.h"

#include <iostream>
#include <fstream>

using namespace espreso;

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	ECFConfiguration ecf(&argc, &argv);
	ECFRedParameters redParameters = ECFReader::read(ecf, ECFReader::configurationFile, ecf.default_args, ecf.variables);

	// add parameters that will be always printed
	redParameters.parameters.push_back(ecf.getParameter(&ecf.input));
	redParameters.parameters.push_back(ecf.getParameter(&ecf.physics));

	auto eachProperties = [&] (ECFObject &object) {
		object.forEachParameters([&] (ECFParameter* parameter) {
			redParameters.parameters.push_back(parameter);
		});
	};

	auto loadStepProperties = [&] (LoadStepConfiguration &loadstep) {
		redParameters.parameters.push_back(loadstep.getParameter(&loadstep.nonlinear_solver));
		redParameters.parameters.push_back(loadstep.getParameter(&loadstep.transient_solver));
		redParameters.parameters.push_back(loadstep.getParameter(&loadstep.duration_time));
		redParameters.parameters.push_back(loadstep.getParameter(&loadstep.type));
		redParameters.parameters.push_back(loadstep.getParameter(&loadstep.mode));

		loadstep.nonlinear_solver.forEachParameters([&] (ECFParameter* parameter) {
			if (
					parameter->data() != &loadstep.nonlinear_solver.adaptive_precision &&
					parameter->data() != &loadstep.nonlinear_solver.c_fact &&
					parameter->data() != &loadstep.nonlinear_solver.r_tol) {

				redParameters.parameters.push_back(parameter);
			}
		});

		eachProperties(loadstep.transient_solver);
	};

	auto materialProperties = [&] (std::map<std::string, MaterialConfiguration> &materials) {
		for (auto mat = materials.begin(); mat != materials.end(); ++mat) {
			mat->second.forEachParameters([&] (ECFParameter* parameter) {
				if (
						parameter->data() != &mat->second.coordinate_system &&
						parameter->data() != &mat->second.physical_model) {
					redParameters.parameters.push_back(parameter);
				}
			});
		}
	};

	auto heatTranferProperties = [&] (HeatTransferConfiguration &configuration) {
		redParameters.parameters.push_back(configuration.getParameter(&configuration.load_steps));
		redParameters.parameters.push_back(configuration.getParameter(&configuration.stabilization));
		redParameters.parameters.push_back(configuration.getParameter(&configuration.sigma));
		for (auto step = configuration.load_steps_settings.begin(); step != configuration.load_steps_settings.end(); ++step) {
			loadStepProperties(step->second);
		}
		materialProperties(configuration.materials);
	};

	auto structuralMechanicsProperties = [&] (StructuralMechanicsConfiguration &configuration) {
		redParameters.parameters.push_back(configuration.getParameter(&configuration.load_steps));
		for (auto step = configuration.load_steps_settings.begin(); step != configuration.load_steps_settings.end(); ++step) {
			loadStepProperties(step->second);
		}
		materialProperties(configuration.materials);
	};

	eachProperties(ecf.workbench);
	eachProperties(ecf.openfoam);
	eachProperties(ecf.esdata);

	heatTranferProperties(ecf.heat_transfer_2d);
	heatTranferProperties(ecf.heat_transfer_3d);
	structuralMechanicsProperties(ecf.structural_mechanics_2d);
	structuralMechanicsProperties(ecf.structural_mechanics_3d);

	redParameters.parameters.push_back(ecf.output.getParameter(&ecf.output.path));
	redParameters.parameters.push_back(ecf.output.getParameter(&ecf.output.format));
	redParameters.parameters.push_back(ecf.output.getParameter(&ecf.output.results_store_frequency));
	redParameters.parameters.push_back(ecf.output.getParameter(&ecf.output.results_nth_stepping));
	redParameters.parameters.push_back(ecf.output.getParameter(&ecf.output.monitors_store_frequency));
	redParameters.parameters.push_back(ecf.output.getParameter(&ecf.output.monitors_nth_stepping));
	redParameters.parameters.push_back(ecf.output.getParameter(&ecf.output.store_results));
	redParameters.parameters.push_back(ecf.output.getParameter(&ecf.output.results_selection));
	eachProperties(ecf.output.results_selection);
	redParameters.parameters.push_back(ecf.output.getParameter(&ecf.output.collected));
	redParameters.parameters.push_back(ecf.output.getParameter(&ecf.output.separate_bodies));
	redParameters.parameters.push_back(ecf.output.getParameter(&ecf.output.separate_materials));
	redParameters.parameters.push_back(ecf.output.getParameter(&ecf.output.domain_shrink_ratio));
	redParameters.parameters.push_back(ecf.output.getParameter(&ecf.output.cluster_shrink_ratio));

	std::ofstream os(ECFReader::configurationFile);
	ECFReader::store(ecf, os, true, false, redParameters);

	ESINFO(OVERVIEW) << "ECF is correct.";

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();

	return 0;
}



