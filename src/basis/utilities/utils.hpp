
#include "../logging/logging.h"
#include <cmath>
#include <algorithm>
#include <cstring>

#include "utils.h"

namespace espreso {

template<typename Ttype>
Ttype Esutils::getEnv(const std::string &name)
{
	Ttype value;
	setFromEnv(value, name);
	return value;
}

template<typename Ttype>
void Esutils::setFromEnv(Ttype &value, const std::string &name)
{
	char *var = getenv(name.c_str());
	if (var != NULL) {
		std::stringstream ss(var);
		ss >> value;
	} else {
		ESINFO(ERROR) << "Set environment variable " << name;
	}
}


template<typename Ttype>
std::vector<Ttype> Esutils::getDistribution(size_t parts, Ttype size)
{
	return getDistribution(parts, (Ttype)0, size);
}

template<typename Ttype>
std::vector<Ttype> Esutils::getDistribution(size_t parts, Ttype start, Ttype end)
{
	if (start > end) {
		ESINFO(ERROR) << "Distribution of interval <" << start << "," << end << "> is not possible.";
	}
	size_t size = end - start;
	std::vector<Ttype> distribution(parts + 1, 0);
	size_t chunkSize = std::ceil(size / (double)parts);
	for (size_t t = 1; t < parts; t++) {
		distribution[t] = t * chunkSize + start;
		if (distribution[t] > size) {
			distribution[t] = end;
		}
	}
	distribution[0] = start;
	distribution[parts] = end;

	return distribution;
}

template<typename Ttype>
Ttype Esutils::sizesToOffsets(std::vector<Ttype> &sizes)
{
	Ttype sum = 0;
	for (size_t i = 0; i < sizes.size(); i++) {
		Ttype tmp = sizes[i];
		sizes[i] = sum;
		sum += tmp;
	}
	return sum;
}

template<typename Ttype>
void Esutils::threadDistributionToFullDistribution(std::vector<std::vector<Ttype> > &distribution)
{
	std::vector<size_t> offsets;
	for (size_t t = 0; t < distribution.size(); t++) {
		offsets.push_back(distribution[t].size() ? distribution[t].back() : 0);
	}
	Esutils::sizesToOffsets(offsets);

	#pragma omp parallel for
	for (size_t t = 0; t < distribution.size(); t++) {
		size_t offset = offsets[t];
		for (size_t i = 0; i < distribution[t].size(); i++) {
			distribution[t][i] += offset;
		}
	}
}

template<typename Ttype>
void Esutils::threadDistributionToFullDistribution(std::vector<Ttype> &data, const std::vector<size_t> &distribution)
{
	size_t threads = distribution.size() - 1;
	std::vector<size_t> offsets(distribution.size());
	for (size_t t = 0; t < threads; t++) {
		if (distribution[t] != distribution[t + 1]) {
			offsets[t] = data[distribution[t + 1] - 1];
		}
	}
	Esutils::sizesToOffsets(offsets);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		size_t offset = offsets[t];
		for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {
			data[i] += offset;
		}
	}
}

template<typename Ttype>
void Esutils::removeDuplicity(std::vector<Ttype> &data, size_t begin)
{
	if (data.size() == begin) {
		return;
	}
	size_t unique = begin;
	for (size_t d = begin + 1; d < data.size(); d++) {
		if (data[unique] != data[d]) {
			data[++unique] = data[d];
		}
	}

	data.resize(unique + 1);
}

template<typename Ttype>
void Esutils::sortAndRemoveDuplicity(std::vector<Ttype> &data, size_t begin)
{
	std::sort(data.begin() + begin, data.end());
	Esutils::removeDuplicity(data, begin);
}

template<typename Ttype>
void Esutils::sortAndRemoveDuplicity(std::vector<std::vector<Ttype> > &data)
{
	for (size_t n = 0; n < data.size(); n++) {
		sortAndRemoveDuplicity(data[n]);
	}
}

template<typename Ttype>
void Esutils::mergeThreadedUniqueData(std::vector<std::vector<Ttype> > &data)
{
	for (size_t t = 1; t < data.size(); t++) {
		data[0].insert(data[0].end(), data[t].begin(), data[t].end());
	}
	sortAndRemoveDuplicity(data[0]);
}

template<typename Ttype>
void Esutils::mergeThreadedUniqueData(std::vector<std::vector<std::vector<Ttype> > > &data)
{
	#pragma omp parallel for
	for (size_t n = 0; n < data[0].size(); n++) {
		for (size_t t = 1; t < data.size(); t++) {
			data[0][n].insert(data[0][n].end(), data[t][n].begin(), data[t][n].end());
		}
		sortAndRemoveDuplicity(data[0][n]);
	}
}

template<typename Ttype>
void Esutils::sortWithInplaceMerge(std::vector<Ttype> &data, const std::vector<size_t> &distribution)
{
	size_t size = distribution.size() - 1;

	#pragma omp parallel for
	for (size_t t = 0; t < size; t++) {
		std::sort(
				data.data() + distribution[t],
				data.data() + distribution[t + 1]);
	}

	size_t align = 1;
	while (align < distribution.size()) align = align << 1;

	std::vector<size_t> _distribution = distribution;
	_distribution.insert(_distribution.end(), align - size, _distribution.back());

	for (size_t i = 2; i <= align; i *= 2) {
		#pragma omp parallel for
		for (size_t t = 0; t < align / i; t++) {
			std::inplace_merge(
					data.data() + _distribution[i * t],
					data.data() + _distribution[i * t + i / 2],
					data.data() + _distribution[i * t + i]);
		}
	}
}

template<typename Ttype>
void Esutils::mergeAppendedData(std::vector<Ttype> &data, const std::vector<size_t> &distribution)
{
	std::vector<size_t> _distribution(distribution.begin() + 1, distribution.end());

	size_t align = 1;
	while (align < _distribution.size()) align = align << 1;

	_distribution.insert(_distribution.end(), align - distribution.size() + 2, _distribution.back());

	for (size_t i = 2; i <= align; i *= 2) {
		#pragma omp parallel for
		for (size_t t = 0; t < align / i; t++) {
			std::inplace_merge(
					data.data() + _distribution[i * t],
					data.data() + _distribution[i * t + i / 2],
					data.data() + _distribution[i * t + i]);
		}
	}

	std::inplace_merge(
			data.data(),
			data.data() + _distribution.front(),
			data.data() + _distribution.back());
}

template<typename Ttype>
typename std::vector<Ttype>::const_iterator Esutils::max_element(const std::vector<Ttype> &elements)
{
	auto max = elements.begin();
	for (size_t i = 1; i < elements.size(); i++) {
		if (*max < elements[i]) {
			max = elements.begin() + i;
		}
	}
	return max;
}

template<>
inline size_t Esutils::packedSize(const std::string &data)
{
	return sizeof(size_t) + data.size();
}

template<typename Ttype>
inline size_t Esutils::packedSize(const Ttype &data)
{
	return sizeof(Ttype);
}

template<>
inline size_t Esutils::packedSize(const std::vector<std::string> &data)
{
	size_t size = sizeof(size_t);
	for (size_t i = 0; i < data.size(); i++) {
		size += packedSize(data[i]);
	}
	return size;
}

template<typename Ttype>
inline size_t Esutils::packedSize(const std::vector<Ttype> &data)
{
	return data.size() * sizeof(Ttype) + sizeof(size_t);
}

template<>
inline void Esutils::pack(const std::string &data, char* &p)
{
	size_t size = data.size();
	memcpy(p, &size, packedSize(size));
	p += packedSize(size);

	memcpy(p, data.data(), data.size());
	p += data.size();
}

template<typename Ttype>
inline void Esutils::pack(const Ttype &data, char* &p)
{
	memcpy(p, &data, packedSize(data));
	p += packedSize(data);
}

template<>
inline void Esutils::pack(const std::vector<std::string> &data, char* &p)
{
	size_t size = data.size();
	memcpy(p, &size, packedSize(size));
	p += packedSize(size);

	for (size_t i = 0; i < data.size(); i++) {
		pack(data[i], p);
	}
}

template<typename Ttype>
inline void Esutils::pack(const std::vector<Ttype> &data, char* &p)
{
	size_t size = data.size();

	memcpy(p, &size, packedSize(size));
	p += packedSize(size);

	if (size) {
		memcpy(p, data.data(), data.size() * sizeof(Ttype));
		p += data.size() * sizeof(Ttype);
	}
}

template<>
inline void Esutils::unpack(std::string &data, const char* &p)
{
	size_t size;
	memcpy(&size, p, packedSize(size));
	p += packedSize(size);
	data = std::string(p, size);
	p += size;
}

template<typename Ttype>
inline void Esutils::unpack(Ttype &data, const char* &p)
{
	memcpy(&data, p, packedSize(data));
	p += packedSize(data);
}

template<>
inline void Esutils::unpack(std::vector<std::string> &data, const char* &p)
{
	size_t size;
	memcpy(&size, p, packedSize(size));
	p += packedSize(size);

	if (size) {
		data.resize(size);
		for (size_t i = 0; i < data.size(); i++) {
			unpack(data[i], p);
		}
	}
}

template<typename Ttype>
inline void Esutils::unpack(std::vector<Ttype> &data, const char* &p)
{
	size_t size;
	memcpy(&size, p, packedSize(size));
	p += packedSize(size);

	if (size) {
		data.resize(size);
		memcpy(data.data(), p, data.size() * sizeof(Ttype));
		p += data.size() * sizeof(Ttype);
	}
}

}
