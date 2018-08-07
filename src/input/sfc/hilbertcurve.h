
#ifndef SRC_INPUT_SFC_HILBERTCURVE_H_
#define SRC_INPUT_SFC_HILBERTCURVE_H_

#include "spacefillingcurve.h"

namespace espreso {

struct HilbertCurve: public SpaceFillingCurve {

	HilbertCurve(size_t dimension, size_t depth, std::vector<Point> &coordinates): SpaceFillingCurve(dimension, depth, coordinates) {};

protected:
	void D1toD2(size_t n, size_t d, size_t &x, size_t &y) const;
	void D1toD3(size_t n, size_t d, size_t &x, size_t &y, size_t &z) const;
	size_t D2toD1(size_t n, size_t x, size_t y) const;
	size_t D3toD1(size_t n, size_t x, size_t y, size_t z) const;

private:
	void rotateD2(size_t n, size_t &x, size_t &y, int rx, int ry) const;
	void rotateD3(size_t n, size_t &x, size_t &y, size_t &z, int rx, int ry, int rz) const;
};

}



#endif /* SRC_INPUT_SFC_HILBERTCURVE_H_ */
