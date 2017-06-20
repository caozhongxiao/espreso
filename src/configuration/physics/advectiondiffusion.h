
#ifndef SRC_CONFIGURATION_PHYSICS_ADVECTIONDIFFUSION_H_
#define SRC_CONFIGURATION_PHYSICS_ADVECTIONDIFFUSION_H_

#include "../material/coordinatesystem.h"
#include "../material/holder.h"
#include "solver.h"
#include "../solver.h"
#include "advectiondiffusionconvection.h"
#include "advectiondiffusionradiation.h"
#include "advectiondiffusionsolver.h"

namespace espreso {

struct AdvectionDiffusionConfiguration: public Configuration {

	enum class STABILIZATION {
		SUPG = 0,
		CAU = 1
	};

	OPTION(STABILIZATION, stabilization, "The type of the stabilization.", STABILIZATION::SUPG, OPTIONS({
		{ "SUPG", STABILIZATION::SUPG, "SUPG stabilization." },
		{ "CAU" , STABILIZATION::CAU , "CAU stabilization." },
	}));
	PARAMETER(double, sigma, "Inconsistent stabilization parameters.", 0);
	PARAMETER(bool, tangent_matrix_correction, "Add derivation matrix to stiffness matrix.", 0);

	PARAMETER(bool, newassembler, "New version of assembler.", 1);

	SUBCONFIG(PhysicsSolver<AdvectionDiffusionNonLinearConvergence>, physics_solver, "Settings of physics solver.");

	SUBMAPTOMAP(size_t, std::string, std::string, heat_flux, "Heat flux", "1", "Heat flux settings for load step '1'", "<REGION>", "<EXPRESSION>");
	SUBMAPTOMAP(size_t, std::string, std::string, heat_flow, "Heat flow", "1", "Heat flow settings for load step '1'", "<REGION>", "<EXPRESSION>");

	SUBMAPTOMAPTOCONFIG(size_t, std::string, AdvectionDiffusionConvection, convection, "Region with convective heat flux",
			"<REGION_NAME>", "Convection for a given region", "1", "Settings for load step '1'");
	SUBMAPTOMAPTOCONFIG(size_t, std::string, AdvectionDiffusionRadiation, diffuse_radiation, "Region with diffuse radiation",
			"<REGION_NAME>", "Diffuse radiation for a given region", "1", "Settings for load step '1'");

	SUBMAP(std::string, std::string, initial_temperature , "Regions initial temperature", "<REGION>", "<EXPRESSION>");

	SUBMAPTOMAP(size_t, std::string, std::string, temperature        , "Temperature"       , "1", "Temperature settings for load step '1'", "<REGION>", "<EXPRESSION>");
	SUBMAPTOMAP(size_t, std::string, std::string, heat_source        , "Heat source"       , "1", "Heat source settings for load step '1'", "<REGION>", "<EXPRESSION>");
	SUBMAPTOMAP(size_t, std::string, std::string, translation_motions, "Translation motion", "1", "Translation motion settings for load step '1'", "<REGION>", "<EXPRESSION>");
	SUBMAPTOMAP(size_t, std::string, std::string, thickness          , "Thickness"         , "1", "Thickness settings for load step '1'", "<REGION>", "<EXPRESSION>");

	SUBMAP(std::string, std::string, material_set, "Assign materials to regions", "<REGION>", "<MATERIAL_NAME>");

	PARAMETER(bool, post_process, "Turn on/off results post processing.", true);

	PARAMETER(bool, bem4i, "Assemble matrices using BEM4I library.", false);
};

}



#endif /* SRC_CONFIGURATION_PHYSICS_ADVECTIONDIFFUSION_H_ */