
#ifndef SRC_PHYSICS_ASSEMBLER_KERNELS_HEATTRANSFER_KERNEL_H_
#define SRC_PHYSICS_ASSEMBLER_KERNELS_HEATTRANSFER_KERNEL_H_

#include "kernel.h"

namespace espreso {

struct ConvectionConfiguration;
struct HeatTransferGlobalSettings;
struct HeatTransferOutputSettings;

struct HeatTransferKernel: public Kernel
{
	HeatTransferKernel(const HeatTransferGlobalSettings &settings, const HeatTransferOutputSettings &output)
	: _settings(settings), _output(output)
	{

	}

	static double convectionHTC(
			const ConvectionConfiguration &convection,
			eslocal csize, double *coordinates, double time, double temp);

protected:
	const HeatTransferGlobalSettings &_settings;
	const HeatTransferOutputSettings &_output;
};

}



#endif /* SRC_PHYSICS_ASSEMBLER_KERNELS_HEATTRANSFER_KERNEL_H_ */
