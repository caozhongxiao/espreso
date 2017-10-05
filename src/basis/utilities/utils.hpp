
#include "../logging/logging.h"
#include <cmath>
#include <algorithm>

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
void Esutils::removeDuplicity(std::vector<Ttype> &data)
{
	if (data.size() == 0) {
		return;
	}
	size_t unique = 0;
	for (size_t d = 1; d < data.size(); d++) {
		if (data[unique] != data[d]) {
			data[++unique] = data[d];
		}
	}

	data.resize(unique + 1);
}

template<typename Ttype>
void Esutils::sortAndRemoveDuplicity(std::vector<Ttype> &data)
{
	std::sort(data.begin(), data.end());
	Esutils::removeDuplicity(data);
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

}
