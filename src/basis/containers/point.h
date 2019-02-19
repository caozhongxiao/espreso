#ifndef POINT3D_H_
#define POINT3D_H_

#include <cmath>

namespace espreso {

class Point {

public:
	Point(): x(0), y(0), z(0) { };
	Point(double x, double y, double z): x(x), y(y), z(z) { };
	Point(const Point &p): x(p.x), y(p.y), z(p.z) { };

	double norm()
	{
		return sqrt(x * x + y * y + z * z);
	}

	void normalize()
	{
		double l = 1.0 / norm();
		x *= l;
		y *= l;
		z *= l;
	}

	double length() const
	{
		return sqrt(x * x + y * y + z * z);
	}

	const Point operator-() const
	{
		return Point(-x, -y, -z);
	}

	void flip()
	{
		x = -x;
		y = -y;
		z = -z;
	}

	const Point operator-(const Point &p) const
	{
		return Point(x - p.x, y - p.y, z - p.z);
	}

	const Point operator+(const Point &p) const
	{
		return Point(x + p.x, y + p.y, z + p.z);
	}

	const Point operator*(double scalar) const
	{
		return Point(x * scalar, y * scalar, z * scalar);
	}

	const Point operator/(double scalar) const
	{
		double xscalar = 1. / scalar;
		return Point(x * xscalar, y * xscalar, z * xscalar);
	}

	double operator*(const Point &p) const
	{
		return x * p.x + y * p.y + z * p.z;
	}

	const Point operator/(const Point &p) const
	{
		return Point(x / p.x, y / p.y, z / p.z);
	}

	static Point cross(const Point &p1, const Point &p2)
	{
		return Point(
				p1.y * p2.z - p1.z * p2.y,
				p1.z * p2.x - p1.x * p2.z,
				p1.x * p2.y - p1.y * p2.x);
	}

	Point &operator*=(double scalar)
	{
		x *= scalar;
		y *= scalar;
		z *= scalar;
		return *this;
	}

	Point &operator+=(const Point &p)
	{
		x += p.x;
		y += p.y;
		z += p.z;
		return *this;
	}

	Point &operator-=(const Point &p)
	{
		x -= p.x;
		y -= p.y;
		z -= p.z;
		return *this;
	}

	Point &operator=(const Point &p)
	{
		x = p.x;
		y = p.y;
		z = p.z;
		return *this;
	}

	Point &operator/=(double scalar)
	{
		x /= scalar;
		y /= scalar;
		z /= scalar;
		return *this;
	}

	double x, y, z;
};

inline bool operator<(const Point &p1, const Point &p2) {
	if (p1.x == p2.x) {
		if (p1.y == p2.y) {
			if (p1.z == p2.z) {
				return false;
			}
			return p1.z < p2.z;
		}
		return p1.y < p2.y;
	}
	return p1.x < p2.x;
}

inline bool operator>(const Point &p1, const Point &p2) {
	return p2 < p1;
}

inline bool operator==(const Point &p1, const Point &p2) {
	return !((p1 < p2) || (p2 < p1));
}

}

#endif /* POINT3D_H_ */
