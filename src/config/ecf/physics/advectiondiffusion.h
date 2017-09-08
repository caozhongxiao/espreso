
#ifndef SRC_CONFIG_ECF_PHYSICS_ADVECTIONDIFFUSION_H_
#define SRC_CONFIG_ECF_PHYSICS_ADVECTIONDIFFUSION_H_

#include "solver.h"
#include "discretization.h"
#include "../material/material.h"

namespace espreso {

struct ConvectionConfiguration: public ECFObject {

	enum class TYPE {
		USER,
		EXTERNAL_NATURAL,
		INTERNAL_NATURAL,
		EXTERNAL_FORCED,
		INTERNAL_FORCED
	};

	enum class VARIANT {
		VERTICAL_WALL,
		INCLINED_WALL,
		HORIZONTAL_CYLINDER,
		SPHERE,
		HORIZONTAL_PLATE_UP,
		HORIZONTAL_PLATE_DOWN,
		AVERAGE_PLATE,
		PARALLEL_PLATES,
		CIRCULAR_TUBE,
		TUBE
	};

	enum class FLUID {
		AIR,
		WATER,
		ENGINE_OIL,
		TRANSFORMER_OIL,
	};

	TYPE type;
	VARIANT variant;
	FLUID fluid;

	std::string heat_transfer_coefficient, external_temperature;
	std::string wall_height, tilt_angle, diameter, plate_length, fluid_velocity, plate_distance, length, absolute_pressure;

	ConvectionConfiguration();
};

struct RadiationConfiguration: public ECFObject {

	std::string emissivity, external_temperature;

	RadiationConfiguration();
};

struct AdvectionDiffusionConfiguration: public ECFObject {

	enum class STABILIZATION {
		SUPG = 0,
		CAU = 1
	};

	STABILIZATION stabilization;
	double sigma;

	AdvectionDiffusionPhysicsSolverConfiguration physics_solver;

	std::map<std::string, std::string> material_set, initial_temperature;
	std::map<size_t, std::map<std::string, std::string> > temperature, heat_source, translation_motions, heat_flux, heat_flow;
	std::map<size_t, std::map<std::string, ConvectionConfiguration> > convection;
	std::map<size_t, std::map<std::string, RadiationConfiguration> > diffuse_radiation;

	bool post_process;

	AdvectionDiffusionConfiguration();
};

struct AdvectionDiffusion2DConfiguration: public AdvectionDiffusionConfiguration {

	std::map<size_t, std::map<std::string, std::string> > thickness;
	std::map<std::string, MaterialConfiguration> materials;

	AdvectionDiffusion2DConfiguration();
};

struct AdvectionDiffusion3DConfiguration: public AdvectionDiffusionConfiguration {

	DISCRETIZATION discretization;
	std::map<std::string, MaterialConfiguration> materials;

	AdvectionDiffusion3DConfiguration();
};

}


#endif /* SRC_CONFIG_ECF_PHYSICS_ADVECTIONDIFFUSION_H_ */
