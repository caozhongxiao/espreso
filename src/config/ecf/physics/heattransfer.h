
#ifndef SRC_CONFIG_ECF_PHYSICS_HEATTRANSFER_H_
#define SRC_CONFIG_ECF_PHYSICS_HEATTRANSFER_H_

#include "physics.h"
#include "physicssolver/loadstep.h"

namespace espreso {

struct ConvectionConfiguration: public ECFDescription {

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
		TUBE,
		QUENCH
	};

	enum class FLUID {
		AIR,
		WATER,
		ENGINE_OIL,
		TRANSFORMER_OIL,
		STEAM
	};

	TYPE type;
	VARIANT variant;
	FLUID fluid;

	ECFExpression heat_transfer_coefficient, external_temperature;
	ECFExpression wall_height, tilt_angle, diameter, plate_length, fluid_velocity, plate_distance, length, experimental_constant, volume_fraction, absolute_pressure;

	ConvectionConfiguration();
};

struct RadiationConfiguration: public ECFDescription {

	ECFExpression emissivity, external_temperature;

	RadiationConfiguration();
};

struct HeatTransferGlobalSettings {

	enum class STABILIZATION {
		SUPG = 0,
		CAU = 1
	};

	STABILIZATION stabilization;
	double sigma;
	bool init_temp_respect_bc, diffusion_split;

	std::map<std::string, ECFExpression> initial_temperature, thickness;
};

struct HeatTransferLoadStepConfiguration: public HeatTransferLoadStepSolverConfiguration {

	std::map<std::string, ECFExpression> temperature, heat_source, heat_flux, heat_flow;
	std::map<std::string, ECFExpressionVector> translation_motions;
	std::map<std::string, ConvectionConfiguration> convection;
	std::map<std::string, RadiationConfiguration> diffuse_radiation;

	HeatTransferLoadStepConfiguration(DIMENSION *D);
};

struct HeatTransferOutputSettings {

	bool temperature, translation_motions, gradient, flux, phase, latent_heat;

	void basic() {
		temperature = translation_motions = true;
		gradient = flux = phase = latent_heat = false;
	}
	void all() {
		temperature = translation_motions = gradient = flux = phase = latent_heat = true;
	}

	HeatTransferOutputSettings() { basic(); }
};

struct HeatTransferConfiguration: public PhysicsConfiguration, public HeatTransferGlobalSettings {

	std::map<size_t, HeatTransferLoadStepConfiguration> load_steps_settings;

	HeatTransferConfiguration(DIMENSION D);
};

}


#endif /* SRC_CONFIG_ECF_PHYSICS_HEATTRANSFER_H_ */
