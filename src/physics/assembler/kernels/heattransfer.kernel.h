
#ifndef SRC_PHYSICS_ASSEMBLER_KERNELS_HEATTRANSFER_KERNEL_H_
#define SRC_PHYSICS_ASSEMBLER_KERNELS_HEATTRANSFER_KERNEL_H_

#include "kernel.h"

namespace espreso {

struct ConvectionConfiguration;

struct HeatTransferKernel: public Kernel
{
	static double convectionHTC(
			const ConvectionConfiguration &convection,
			esint csize, double *coordinates, double time, double temp);
};

}



#endif /* SRC_PHYSICS_ASSEMBLER_KERNELS_HEATTRANSFER_KERNEL_H_ */
