
#ifndef SRC_CONFIG_ADVECTIONDIFFUSION2D_H_
#define SRC_CONFIG_ADVECTIONDIFFUSION2D_H_

#include "solver.h"
#include "coordinatesystem.h"
#include "advectiondiffusionconvection.h"
#include "advectiondiffusionsolver.h"
#include "materialholder.h"

namespace espreso {

struct AdvectionDiffusion2DMaterial: public Configuration {

	PARAMETER(MaterialParam<MATERIAL_PARAMETER::DENSITY>                , density, "Density", {"0"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::HEAT_CAPACITY>          , Cp     , "Termal capacity."       , {"0"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XX>, KXX    , "Termal conductivity XX.", {"1"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_YY>, KYY    , "Termal conductivity YY.", {"1"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XY>, KXY    , "Termal conductivity XY.", {"1"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_YX>, KYX    , "Termal conductivity YX.", {"1"});

	OPTION(MATERIAL_MODEL, model, "Material model", MATERIAL_MODEL::ISOTROPIC, OPTIONS({
		{ "ISOTROPIC"  , MATERIAL_MODEL::ISOTROPIC  , "Isotropic." },
		{ "DIAGONAL"   , MATERIAL_MODEL::DIAGONAL   , "Diagonal." },
		{ "SYMMETRIC"  , MATERIAL_MODEL::SYMMETRIC  , "Symmetric." },
		{ "ANISOTROPIC", MATERIAL_MODEL::ANISOTROPIC, "Anisotropic." }
	}));

	SUBCONFIG(CoordinateSystem, coordinate_system, "Element coordinate system.");
};

struct AdvectionDiffusion2DConfiguration: public Configuration {

	enum class STABILIZATION {
		SUPG = 0,
		CAU = 1
	};

	OPTION(STABILIZATION, stabilization, "The type of the stabilization.", STABILIZATION::SUPG, OPTIONS({
		{ "SUPG", STABILIZATION::SUPG, "SUPG stabilization." },
		{ "CAU" , STABILIZATION::CAU , "CAU stabilization." },
	}));
	PARAMETER(double, sigma, "Inconsistent stabilization parameters.", 0);

	SUBCONFIG(AdvectionDiffusionSolver, physics_solver, "Settings of physics solver.");

	OPTION(SOLVER_LIBRARY, solver_library, "Linear solver used for computing a system.", SOLVER_LIBRARY::ESPRESO, OPTIONS({
		{ "ESPRESO", SOLVER_LIBRARY::ESPRESO, "ESPRESO solver [FETI methods]" },
		{ "HYPRE"  , SOLVER_LIBRARY::HYPRE  , "Hypre solver [multigrid methods]" },
	}));

	SUBCONFIG(ESPRESOSolver, espreso, "Internal FETI solver options.");
	SUBCONFIG(HypreSolver  , hypre  , "Multigrid solver setting.");

	SUBVECTORMAP(std::string, std::string, heat_flux, "Heat flux", "<LOAD_STEPS>", "<REGION>", "<EXPRESSION>");
	SUBVECTORMAP(std::string, std::string, heat_flow, "Heat flow", "<LOAD_STEPS>", "<REGION>", "<EXPRESSION>");

	SUBMAPTOMAPTOCONFIG(size_t, std::string, AdvectionDiffusionConvection, convection, "Region with convective heat flux");

	SUBMAP(std::string, std::string, initial_temperature , "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");

	SUBVECTORMAP(std::string, std::string, temperature        , "Temperature"       , "<LOAD_STEPS>", "<REGION>", "<EXPRESSION>");
	SUBVECTORMAP(std::string, std::string, heat_source        , "Heat source"       , "<LOAD_STEPS>", "<REGION>", "<EXPRESSION>");
	SUBVECTORMAP(std::string, std::string, translation_motions, "Translation motion", "<LOAD_STEPS>", "<REGION>", "<EXPRESSION>");
	SUBVECTORMAP(std::string, std::string, thickness          , "Thicness"          , "<LOAD_STEPS>", "<REGION>", "<EXPRESSION>");

	SUBMAPTOCONFIG(std::string, AdvectionDiffusion2DMaterial, materials, "Material description.");
	SUBMAP(std::string, std::string, material_set, "Assign materials to regions", "<REGION>", "<MATERIAL_NAME>");

	std::map<std::string, std::map<std::string, AdvectionDiffusion2DMaterial*> > xxx = mapToMapToConfiguration<std::string, std::string, AdvectionDiffusion2DMaterial>::create("xxx", "DESC", this);

	PARAMETER(bool, post_process, "Turn on/off results post processing.", true);

};

}



#endif /* SRC_CONFIG_ADVECTIONDIFFUSION2D_H_ */
