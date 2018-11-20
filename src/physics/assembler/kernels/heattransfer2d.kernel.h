
#ifndef SRC_PHYSICS_ASSEMBLER_KERNELS_HEATTRANSFER2D_KERNEL_H_
#define SRC_PHYSICS_ASSEMBLER_KERNELS_HEATTRANSFER2D_KERNEL_H_

#include "kernel.h"

namespace espreso {

struct DenseMatrix;
struct Element;
struct Step;
struct HeatTransferGlobalSettings;
struct MaterialBaseConfiguration;
struct MaterialConfiguration;
enum Matrices: int;

struct HeatTransfer2DKernel: public Kernel
{
	struct Iterator {
		Element *element;

		double *temperature;
		double *coordinates;
		double *gradient;
		double *motion;
		double *heat;
		double *thickness;

		const MaterialConfiguration *material;
	};

	HeatTransfer2DKernel(const HeatTransferGlobalSettings &settings);

	void processElement(Matrices matrices, const Iterator &iterator, const Step &step, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const;
//	void processFace(eslocal domain, const BoundaryRegionStore *region, Matrices matrices, eslocal findex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const;
//	void processEdge(eslocal domain, const BoundaryRegionStore *region, Matrices matrices, eslocal eindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const;
//	void processNode(eslocal domain, const BoundaryRegionStore *region, Matrices matrices, eslocal nindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const;
//	void processSolution();

protected:
	void assembleMaterialMatrix(eslocal node, double *coordinates, const MaterialBaseConfiguration *mat, double phase, double time, double temp, DenseMatrix &K, DenseMatrix &CD, bool tangentCorrection) const;
//	void postProcessElement(eslocal domain, eslocal eindex);

	const HeatTransferGlobalSettings &_settings;
};

}



#endif /* SRC_PHYSICS_ASSEMBLER_KERNELS_HEATTRANSFER2D_KERNEL_H_ */
