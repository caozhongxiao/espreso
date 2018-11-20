
#ifndef SRC_PHYSICS_ASSEMBLER_KERNELS_KERNEL_H_
#define SRC_PHYSICS_ASSEMBLER_KERNELS_KERNEL_H_

#include <cstddef>

namespace espreso {

struct Kernel {

	void smoothstep(double &smoothStep, double &derivation, double edge0, double edge1, double value, size_t order) const;

	double determinant2x2(double *values) const
	{
		return values[0] * values[3] - values[1] * values[2];
	}

	void inverse2x2(const double *m, double *inv, double det) const
	{
		double detJx = 1 / det;
		inv[0] =   detJx * m[3];
		inv[1] = - detJx * m[1];
		inv[2] = - detJx * m[2];
		inv[3] =   detJx * m[0];
	}

	double determinant3x3(double *values) const
	{
		return
			+ values[0] * values[4] * values[8]
			+ values[1] * values[5] * values[6]
			+ values[2] * values[3] * values[7]
			- values[2] * values[4] * values[6]
			- values[1] * values[3] * values[8]
			- values[0] * values[5] * values[7];
	}

	void inverse3x3(const double *m, double *inv, double det) const
	{
		double detJx = 1 / det;
		inv[0] = detJx * ( m[8] * m[4] - m[7] * m[5]);
		inv[1] = detJx * (-m[8] * m[1] + m[7] * m[2]);
		inv[2] = detJx * ( m[5] * m[1] - m[4] * m[2]);
		inv[3] = detJx * (-m[8] * m[3] + m[6] * m[5]);
		inv[4] = detJx * ( m[8] * m[0] - m[6] * m[2]);
		inv[5] = detJx * (-m[5] * m[0] + m[3] * m[2]);
		inv[6] = detJx * ( m[7] * m[3] - m[6] * m[4]);
		inv[7] = detJx * (-m[7] * m[0] + m[6] * m[1]);
		inv[8] = detJx * ( m[4] * m[0] - m[3] * m[1]);
	}
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_KERNELS_KERNEL_H_ */
