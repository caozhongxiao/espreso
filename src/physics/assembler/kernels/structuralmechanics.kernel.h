
#ifndef SRC_PHYSICS_ASSEMBLER_KERNELS_STRUCTURALMECHANICS_KERNEL_H_
#define SRC_PHYSICS_ASSEMBLER_KERNELS_STRUCTURALMECHANICS_KERNEL_H_

#include "kernel.h"

namespace espreso {

struct StructuralMechanicsGlobalSettings;
struct StructuralMechanicsOutputSettings;

struct StructuralMechanicsKernel: public Kernel
{
	StructuralMechanicsKernel(const StructuralMechanicsGlobalSettings &settings, const StructuralMechanicsOutputSettings &output)
	: _settings(settings), _output(output)
	{

	}

protected:
	const StructuralMechanicsGlobalSettings &_settings;
	const StructuralMechanicsOutputSettings &_output;
};

}



#endif /* SRC_PHYSICS_ASSEMBLER_KERNELS_STRUCTURALMECHANICS_KERNEL_H_ */
