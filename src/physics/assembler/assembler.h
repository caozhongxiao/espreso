
#ifndef SRC_PHYSICS_ASSEMBLER_ASSEMBLER_H_
#define SRC_PHYSICS_ASSEMBLER_ASSEMBLER_H_

#include "composer/global/uniformnodescomposer.h"

#include "controllers/heattransfer2d.controller.h"
#include "controllers/heattransfer3d.controller.h"
#include "controllers/structuralmechanics2d.controller.h"
#include "controllers/structuralmechanics3d.controller.h"

#include "../../config/ecf/physics/heattransfer.h"
#include "../../config/ecf/physics/structuralmechanics.h"
#include "composer/feti/uniformnodefeticomposer.h"

namespace espreso {

struct Assembler {
	void run();
};

template <typename TComposer, typename TController>
struct GlobalAssembler: public Assembler, public TComposer, public TController {

	GlobalAssembler()
	: TComposer(),
	  TController() {}
};

struct DomainsHeatTransfer2D: public HeatTransfer2DControler, public UniformNodeFETIComposer {

	DomainsHeatTransfer2D(
			Mesh &mesh, DataHolder &instance,
			const HeatTransferGlobalSettings &gSettings, const HeatTransferLoadStepConfiguration &sSettings, const HeatTransferOutputSettings &oSettings)

	: HeatTransfer2DControler(mesh, gSettings, sSettings, oSettings),
	  UniformNodeFETIComposer(mesh, instance, *this, 1, sSettings.feti.redundant_lagrange, sSettings.feti.scaling) {}
};

struct GlobalHeatTransfer2D: public HeatTransfer2DControler, public UniformNodesComposer {

	GlobalHeatTransfer2D(
			Mesh &mesh, DataHolder &instance,
			const HeatTransferGlobalSettings &gSettings, const HeatTransferLoadStepConfiguration &sSettings, const HeatTransferOutputSettings &oSettings)

	: HeatTransfer2DControler(mesh, gSettings, sSettings, oSettings), UniformNodesComposer(mesh, instance, *this, 1) {}
};


struct DomainsHeatTransfer3D: public HeatTransfer3DControler, public UniformNodeFETIComposer {

	DomainsHeatTransfer3D(
			Mesh &mesh, DataHolder &instance,
			const HeatTransferGlobalSettings &gSettings, const HeatTransferLoadStepConfiguration &sSettings, const HeatTransferOutputSettings &oSettings)

	: HeatTransfer3DControler(mesh, gSettings, sSettings, oSettings),
	  UniformNodeFETIComposer(mesh, instance, *this, 1, sSettings.feti.redundant_lagrange, sSettings.feti.scaling) {}
};

struct GlobalHeatTransfer3D: public HeatTransfer3DControler, public UniformNodesComposer {

	GlobalHeatTransfer3D(
			Mesh &mesh, DataHolder &instance,
			const HeatTransferGlobalSettings &gSettings, const HeatTransferLoadStepConfiguration &sSettings, const HeatTransferOutputSettings &oSettings)

	: HeatTransfer3DControler(mesh, gSettings, sSettings, oSettings), UniformNodesComposer(mesh, instance, *this, 1) {}
};

struct DomainsStructuralMechanics2D: public StructuralMechanics2DControler, public UniformNodeFETIComposer {

	DomainsStructuralMechanics2D(
			Mesh &mesh, DataHolder &instance,
			const StructuralMechanicsGlobalSettings &gSettings, const StructuralMechanicsLoadStepConfiguration &sSettings, const StructuralMechanicsOutputSettings &oSettings)

	: StructuralMechanics2DControler(mesh, gSettings, sSettings, oSettings),
	  UniformNodeFETIComposer(mesh, instance, *this, 1, sSettings.feti.redundant_lagrange, sSettings.feti.scaling) {}
};

struct GlobalStructuralMechanics2D: public StructuralMechanics2DControler, public UniformNodesComposer {

	GlobalStructuralMechanics2D(
			Mesh &mesh, DataHolder &instance,
			const StructuralMechanicsGlobalSettings &gSettings, const StructuralMechanicsLoadStepConfiguration &sSettings, const StructuralMechanicsOutputSettings &oSettings)

	: StructuralMechanics2DControler(mesh, gSettings, sSettings, oSettings), UniformNodesComposer(mesh, instance, *this, 1) {}
};


struct DomainsStructuralMechanics3D: public StructuralMechanics3DControler, public UniformNodeFETIComposer {

	DomainsStructuralMechanics3D(
			Mesh &mesh, DataHolder &instance,
			const StructuralMechanicsGlobalSettings &gSettings, const StructuralMechanicsLoadStepConfiguration &sSettings, const StructuralMechanicsOutputSettings &oSettings)

	: StructuralMechanics3DControler(mesh, gSettings, sSettings, oSettings),
	  UniformNodeFETIComposer(mesh, instance, *this, 1, sSettings.feti.redundant_lagrange, sSettings.feti.scaling) {}
};

struct GlobalStructuralMechanics3D: public StructuralMechanics3DControler, public UniformNodesComposer {

	GlobalStructuralMechanics3D(
			Mesh &mesh, DataHolder &instance,
			const StructuralMechanicsGlobalSettings &gSettings, const StructuralMechanicsLoadStepConfiguration &sSettings, const StructuralMechanicsOutputSettings &oSettings)

	: StructuralMechanics3DControler(mesh, gSettings, sSettings, oSettings), UniformNodesComposer(mesh, instance, *this, 1) {}
};

}



#endif /* SRC_PHYSICS_ASSEMBLER_ASSEMBLER_H_ */
