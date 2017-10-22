
#ifndef SRC_BASIS_CONTAINERS_TARRAY_H_
#define SRC_BASIS_CONTAINERS_TARRAY_H_

#include <cmath>
#include <vector>
#include <cstring>

#include <omp.h>

namespace espreso {


/// Threaded array
template <typename TType>
class tarray {

public:
	static std::vector<size_t> distribute(size_t threads, size_t size);

	tarray(size_t threads, size_t size);
	tarray(const std::vector<std::vector<TType> > &data);
	tarray(size_t thread, size_t threads, const std::vector<TType> &data);

	tarray(const tarray<TType> &other);
	tarray(tarray<TType> &&other);
	tarray<TType>& operator=(const tarray<TType> &other);
	tarray<TType>& operator=(tarray<TType> &&other);

	size_t                     size()         const { return _size; }
	size_t                     threads()      const { return _distribution.size() - 1; }
	const std::vector<size_t>& distribution() const { return _distribution; }

	TType& operator[] (size_t n) { return _data[n]; }
	TType* data()                { return _data; }
	TType& back()                { return _data[_size - 1]; }
	TType& front()               { return _data[0]; }
	TType* begin()               { return _data; }
	TType* end()                 { return _data + _size; }

	const TType& operator[] (size_t n) const { return _data[n]; }
	const TType* data()                const { return _data; }
	const TType& back()                const { return _data[_size - 1]; }
	const TType& front()               const { return _data[0]; }
	const TType* begin()               const { return _data; }
	const TType* end()                 const { return _data + _size; }
	const TType* cbegin()              const { return _data; }
	const TType* cend()                const { return _data + _size; }

	~tarray();

private:
	size_t _size;
	TType *_data;

	std::vector<size_t> _distribution;
};

}

#include "tarray.hpp"

#endif /* SRC_BASIS_CONTAINERS_TARRAY_H_ */
