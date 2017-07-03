
#include "assembler.h"
#include "old_physics/elasticity2d/assembler.h"
#include "old_physics/elasticity3d/assembler.h"
#include "old_physics/linear/advectiondiffusion2d/assembler.h"
#include "old_physics/linear/advectiondiffusion3d/assembler.h"
#include "old_physics/linear/laplacesteklovpoincare/assembler.h"

#include <numeric>

#include "../configuration/globalconfiguration.h"
#include "instance/dynamics/instance.h"
#include "instance/hypre/instance.h"
#include "instance/linear/instance.h"
#include "instance/nonlinear/ssnm/instance.h"

using namespace espreso;

class Mesh;

void Assembler::compose(const GlobalConfiguration &configuration, OldInstance* &instance, Mesh &mesh)
{
	switch (configuration.physics) {
	case PHYSICS::STRUCTURAL_MECHANICS_2D:
	case PHYSICS::LINEAR_ELASTICITY_2D:
		switch (configuration.linear_elasticity_2D.physics_solver.load_steps_settings.at(1)->solver_library) {
		case SOLVER_LIBRARY::ESPRESO:
			instance = new LinearInstance<Elasticity2D, LinearElasticity2DConfiguration>(configuration.linear_elasticity_2D, configuration.output, mesh);
			break;
		case SOLVER_LIBRARY::HYPRE:
			instance = new HypreInstance<Elasticity2D, LinearElasticity2DConfiguration>(configuration.linear_elasticity_2D, configuration.output, mesh);
			break;
		};
		break;
	case PHYSICS::STRUCTURAL_MECHANICS_3D:
	case PHYSICS::LINEAR_ELASTICITY_3D:
		switch (configuration.linear_elasticity_3D.physics_solver.load_steps_settings.at(1)->solver_library) {
		case SOLVER_LIBRARY::ESPRESO:
			instance = new LinearInstance<Elasticity3D, LinearElasticity3DConfiguration>(configuration.linear_elasticity_3D, configuration.output, mesh);
			break;
		case SOLVER_LIBRARY::HYPRE:
			instance = new HypreInstance<Elasticity3D, LinearElasticity3DConfiguration>(configuration.linear_elasticity_3D, configuration.output, mesh);
			break;
		};
		break;
	case PHYSICS::ADVECTION_DIFFUSION_2D:
		switch (configuration.advection_diffusion_2D.physics_solver.load_steps_settings.at(1)->solver_library) {
		case SOLVER_LIBRARY::ESPRESO:
			instance = new LinearInstance<AdvectionDiffusion2D, AdvectionDiffusion2DConfiguration>(configuration.advection_diffusion_2D, configuration.output, mesh);
			break;
		case SOLVER_LIBRARY::HYPRE:
			instance = new HypreInstance<AdvectionDiffusion2D, AdvectionDiffusion2DConfiguration>(configuration.advection_diffusion_2D, configuration.output, mesh);
			break;
		};
		break;
	case PHYSICS::ADVECTION_DIFFUSION_3D:
		switch (configuration.advection_diffusion_3D.physics_solver.load_steps_settings.at(1)->solver_library) {
		case SOLVER_LIBRARY::ESPRESO:
			if (configuration.advection_diffusion_3D.bem4i) {
				instance = new LinearInstance<LaplaceSteklovPoincare, AdvectionDiffusion3DConfiguration>(configuration.advection_diffusion_3D, configuration.output, mesh);
			} else {
				instance = new LinearInstance<AdvectionDiffusion3D, AdvectionDiffusion3DConfiguration>(configuration.advection_diffusion_3D, configuration.output, mesh);
			}
			break;
		case SOLVER_LIBRARY::HYPRE:
			instance = new HypreInstance<AdvectionDiffusion3D, AdvectionDiffusion3DConfiguration>(configuration.advection_diffusion_3D, configuration.output, mesh);
			break;
		};
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Not implemented physics";
	}

}



