
#ifndef SRC_BASIS_CONTAINERS_TARRAY_HPP_
#define SRC_BASIS_CONTAINERS_TARRAY_HPP_

#include "../containers/tarray.h"

namespace espreso {

template <typename TType>
std::vector<size_t> tarray<TType>::distribute(size_t threads, size_t size)
{
	std::vector<size_t> distribution(threads + 1);

	size_t chunkSize = std::ceil(size / (double)threads);
	for (size_t t = 1; t < threads; t++) {
		distribution[t] = t * chunkSize;
		if (distribution[t] > size) {
			distribution[t] = size;
		}
	}
	distribution[threads] = size;

	return distribution;
}


template <typename TType>
tarray<TType>::tarray(const tarray<TType> &other)
: _size(other._size), _data(other._data), _distribution(other._distribution)
{
	if (_size) {
		_data = new TType[_size];
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			memcpy(_data + _distribution[t], other._data + _distribution[t], (_distribution[t + 1] - _distribution[t]) * sizeof(TType));
		}
	}
}

template <typename TType>
tarray<TType>::tarray(tarray<TType> &&other)
: _size(std::move(other._size)), _data(std::move(other._data)), _distribution(std::move(other._distribution))
{
	other._size = 0;
	other._data = NULL;
}

template <typename TType>
tarray<TType>& tarray<TType>::operator=(const tarray<TType> &other)
{
	if (this != &other) {
		_size = other._size;
		_distribution = other._distribution;

		if (_data != NULL) {
			delete[] _data;
		}
		_data = new TType[_size];
		#pragma omp parallel for
		for (size_t t = 0; t < other.threads(); t++) {
			memcpy(_data + _distribution[t], other._data + _distribution[t], (_distribution[t + 1] - _distribution[t]) * sizeof(TType));
		}
	}
	return *this;
}

template <typename TType>
tarray<TType>& tarray<TType>::operator=(tarray<TType> &&other)
{
	if (this != &other) {
		_size = std::move(other._size);
		_data = std::move(other._data);
		_distribution = std::move(other._distribution);

		other._size = 0;
		other._data = NULL;
	}
	return *this;
}

template <typename TType>
tarray<TType>::tarray(size_t threads, size_t size)
:  _size(size), _data(NULL)
{
	_distribution = distribute(threads, size);

	if (_size) {
		_data = new TType[_size];
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (size_t i = _distribution[t]; i < _distribution[t + 1]; i++) {
				_data[i] = TType{};
			}
		}
	}
}

template <typename TType>
tarray<TType>::tarray(const std::vector<std::vector<TType> > &data)
: _size(0), _data(NULL)
{
	_distribution = std::vector<size_t>(data.size() + 1, 0);
	for (size_t t = 0; t < data.size(); t++) {
		_size += data[t].size();
		_distribution[t + 1] = _size;
	}

	if (_size) {
		_data = new TType[_size];
		#pragma omp parallel for
		for (size_t t = 0; t < threads(); t++) {
			memcpy(_data + _distribution[t], data[t].data(), data[t].size() * sizeof(TType));
		}
	}
}

template <typename TType>
tarray<TType>::~tarray()
{
	if (_data) {
		delete[] _data;
	}
}

}

#endif /* SRC_BASIS_CONTAINERS_TARRAY_HPP_ */
