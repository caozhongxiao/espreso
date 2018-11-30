
#ifndef SRC_PHYSICS_ASSEMBLER_ASSEMBLER_H_
#define SRC_PHYSICS_ASSEMBLER_ASSEMBLER_H_

#include "composer/domains/uniformnodedomainscomposer.h"
#include "composer/global/uniformnodescomposer.h"

#include "controllers/heattransfer2d.controller.h"

namespace espreso {

struct LoadStepConfiguration;
struct HeatTransferLoadStepConfiguration;

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


}



#endif /* SRC_PHYSICS_ASSEMBLER_ASSEMBLER_H_ */
