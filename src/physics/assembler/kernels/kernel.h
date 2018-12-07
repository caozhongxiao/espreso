
#ifndef SRC_PHYSICS_ASSEMBLER_KERNELS_KERNEL_H_
#define SRC_PHYSICS_ASSEMBLER_KERNELS_KERNEL_H_

#include <cstddef>
#include <cstring>

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

	// source dX, dY

	// target::
	// dX   0
	//  0  dY
	// dY  dX
	void distribute3x2(double *target, double *source, size_t rows, size_t columns) const
	{
		memcpy(target                               , source          , sizeof(double) * columns);
		memcpy(target + 2 * rows * columns + columns, source          , sizeof(double) * columns);

		memcpy(target + 1 * rows * columns + columns, source + columns, sizeof(double) * columns);
		memcpy(target + 2 * rows * columns          , source + columns, sizeof(double) * columns);
	}

	// source dX, dY

	// target::
	// dX   0
	//  0  dY
	//  0  0
	// dY  dX
	void distribute4x2(double *target, double *source, size_t rows, size_t columns) const
	{
		memcpy(target                               , source          , sizeof(double) * columns);
		memcpy(target + 3 * rows * columns + columns, source          , sizeof(double) * columns);

		memcpy(target + 1 * rows * columns + columns, source + columns, sizeof(double) * columns);
		memcpy(target + 3 * rows * columns          , source + columns, sizeof(double) * columns);
	}

	// source dX, dY, dZ

	// target::
	// dX   0   0
	//  0  dY   0
	//  0   0  dZ
	// dY  dX   0
	//  0  dZ  dY
	// dZ   0  dX
	void distribute6x3(double *target, const double *source, size_t rows, size_t columns) const
	{
		const double *dNDx = source;
		const double *dNDy = source + columns;
		const double *dNDz = source + 2 * columns;

		memcpy(target                                   , dNDx, sizeof(double) * columns);
		memcpy(target + 3 * rows * columns +     columns, dNDx, sizeof(double) * columns);
		memcpy(target + 5 * rows * columns + 2 * columns, dNDx, sizeof(double) * columns);

		memcpy(target + 1 * rows * columns +     columns, dNDy, sizeof(double) * columns);
		memcpy(target + 3 * rows * columns              , dNDy, sizeof(double) * columns);
		memcpy(target + 4 * rows * columns + 2 * columns, dNDy, sizeof(double) * columns);

		memcpy(target + 2 * rows * columns + 2 * columns, dNDz, sizeof(double) * columns);
		memcpy(target + 4 * rows * columns +     columns, dNDz, sizeof(double) * columns);
		memcpy(target + 5 * rows * columns              , dNDz, sizeof(double) * columns);
	}
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_KERNELS_KERNEL_H_ */
