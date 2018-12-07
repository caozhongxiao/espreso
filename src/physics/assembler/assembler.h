
#ifndef SRC_PHYSICS_ASSEMBLER_ASSEMBLER_H_
#define SRC_PHYSICS_ASSEMBLER_ASSEMBLER_H_

#include "composer/domains/uniformnodedomainscomposer.h"
#include "composer/global/uniformnodescomposer.h"

#include "controllers/heattransfer2d.controller.h"
#include "controllers/heattransfer3d.controller.h"
#include "controllers/structuralmechanics2d.controller.h"
#include "controllers/structuralmechanics3d.controller.h"

namespace espreso {

struct LoadStepConfiguration;
struct HeatTransferLoadStepConfiguration;
struct StructuralMechanicsLoadStepConfiguration;

struct DomainsHeatTransfer2D: public HeatTransfer2DControler, public UniformNodeDomainsComposer {

	DomainsHeatTransfer2D(
			Mesh &mesh, Instance &instance, Step &step,
			const HeatTransferGlobalSettings &gSettings, const HeatTransferLoadStepConfiguration &sSettings, const HeatTransferOutputSettings &oSettings)

	: HeatTransfer2DControler(mesh, step, gSettings, sSettings, oSettings),
	  UniformNodeDomainsComposer(mesh, step, instance, *this, 1, sSettings.feti.redundant_lagrange, sSettings.feti.scaling) {}
};

struct GlobalHeatTransfer2D: public HeatTransfer2DControler, public UniformNodesComposer {

	GlobalHeatTransfer2D(
			Mesh &mesh, Instance &instance, Step &step,
			const HeatTransferGlobalSettings &gSettings, const HeatTransferLoadStepConfiguration &sSettings, const HeatTransferOutputSettings &oSettings)

	: HeatTransfer2DControler(mesh, step, gSettings, sSettings, oSettings), UniformNodesComposer(mesh, step, instance, *this, 1) {}
};


struct DomainsHeatTransfer3D: public HeatTransfer3DControler, public UniformNodeDomainsComposer {

	DomainsHeatTransfer3D(
			Mesh &mesh, Instance &instance, Step &step,
			const HeatTransferGlobalSettings &gSettings, const HeatTransferLoadStepConfiguration &sSettings, const HeatTransferOutputSettings &oSettings)

	: HeatTransfer3DControler(mesh, step, gSettings, sSettings, oSettings),
	  UniformNodeDomainsComposer(mesh, step, instance, *this, 1, sSettings.feti.redundant_lagrange, sSettings.feti.scaling) {}
};

struct GlobalHeatTransfer3D: public HeatTransfer3DControler, public UniformNodesComposer {

	GlobalHeatTransfer3D(
			Mesh &mesh, Instance &instance, Step &step,
			const HeatTransferGlobalSettings &gSettings, const HeatTransferLoadStepConfiguration &sSettings, const HeatTransferOutputSettings &oSettings)

	: HeatTransfer3DControler(mesh, step, gSettings, sSettings, oSettings), UniformNodesComposer(mesh, step, instance, *this, 1) {}
};

struct DomainsStructuralMechanics2D: public StructuralMechanics2DControler, public UniformNodeDomainsComposer {

	DomainsStructuralMechanics2D(
			Mesh &mesh, Instance &instance, Step &step,
			const StructuralMechanicsGlobalSettings &gSettings, const StructuralMechanicsLoadStepConfiguration &sSettings, const StructuralMechanicsOutputSettings &oSettings)

	: StructuralMechanics2DControler(mesh, step, gSettings, sSettings, oSettings),
	  UniformNodeDomainsComposer(mesh, step, instance, *this, 1, sSettings.feti.redundant_lagrange, sSettings.feti.scaling) {}
};

struct GlobalStructuralMechanics2D: public StructuralMechanics2DControler, public UniformNodesComposer {

	GlobalStructuralMechanics2D(
			Mesh &mesh, Instance &instance, Step &step,
			const StructuralMechanicsGlobalSettings &gSettings, const StructuralMechanicsLoadStepConfiguration &sSettings, const StructuralMechanicsOutputSettings &oSettings)

	: StructuralMechanics2DControler(mesh, step, gSettings, sSettings, oSettings), UniformNodesComposer(mesh, step, instance, *this, 1) {}
};


struct DomainsStructuralMechanics3D: public StructuralMechanics3DControler, public UniformNodeDomainsComposer {

	DomainsStructuralMechanics3D(
			Mesh &mesh, Instance &instance, Step &step,
			const StructuralMechanicsGlobalSettings &gSettings, const StructuralMechanicsLoadStepConfiguration &sSettings, const StructuralMechanicsOutputSettings &oSettings)

	: StructuralMechanics3DControler(mesh, step, gSettings, sSettings, oSettings),
	  UniformNodeDomainsComposer(mesh, step, instance, *this, 1, sSettings.feti.redundant_lagrange, sSettings.feti.scaling) {}
};

struct GlobalStructuralMechanics3D: public StructuralMechanics3DControler, public UniformNodesComposer {

	GlobalStructuralMechanics3D(
			Mesh &mesh, Instance &instance, Step &step,
			const StructuralMechanicsGlobalSettings &gSettings, const StructuralMechanicsLoadStepConfiguration &sSettings, const StructuralMechanicsOutputSettings &oSettings)

	: StructuralMechanics3DControler(mesh, step, gSettings, sSettings, oSettings), UniformNodesComposer(mesh, step, instance, *this, 1) {}
};

}



#endif /* SRC_PHYSICS_ASSEMBLER_ASSEMBLER_H_ */
