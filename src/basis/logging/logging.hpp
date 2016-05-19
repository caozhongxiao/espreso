
#include "logging.h"

namespace espreso {

template<typename Tvalue>
static void gather(const Tvalue &value, Tvalue &min, Tvalue &max, Tvalue &total)
{
	typename std::vector<Tvalue> values(config::env::MPIsize);

	MPI_Gather(&value, sizeof(Tvalue), MPI_BYTE, values.data(), sizeof(Tvalue), MPI_BYTE, 0, MPI_COMM_WORLD);

	min = max = value;
	total = 0;
	for (size_t i = 0; i < values.size(); i++) {
		total += values[i];
		min = std::min(min, values[i]);
		max = std::max(max, values[i]);
	}
}

template<typename Tvalue>
static void gather(const std::pair<Tvalue, Tvalue> &value, Tvalue &min, Tvalue &max)
{
	typename std::vector<std::pair<Tvalue, Tvalue> > values(config::env::MPIsize);

	MPI_Gather(&value, 2 * sizeof(Tvalue), MPI_BYTE, values.data(), 2 * sizeof(Tvalue), MPI_BYTE, 0, MPI_COMM_WORLD);

	min = value.first;
	max = value.second;
	for (size_t i = 0; i < values.size(); i++) {
		min = std::min(min, values[i].first);
		max = std::max(max, values[i].second);
	}
}


template<typename Tvalue>
std::string Info::sumValue(const Tvalue &value)
{
	Tvalue min, max, total;
	gather(value, min, max, total);

	std::stringstream ss;
	ss << total << ", average: " << (double)total / config::env::MPIsize << " (from " << min << " to " << max << ")";
	return ss.str();
}


template<typename Tvalue>
std::string Info::averageValue(const Tvalue &value)
{
	Tvalue min, max, total;
	gather(value, min, max, total);

	std::stringstream ss;
	ss << "average: " << (double)total / config::env::MPIsize << " (from " << min << " to " << max << ")";
	return ss.str();
}

template<typename Tvalue>
std::string Info::averageValues(const std::vector<Tvalue> &values)
{
	std::pair<Tvalue, Tvalue> minmax(values[0], values[0]);
	for (size_t i = 1; i < values.size(); i++) {
		minmax.first = std::min(minmax.first, values[i]);
		minmax.second = std::max(minmax.second, values[i]);
	}

	Tvalue min, max;
	gather(minmax, min, max);

	std::stringstream ss;
	ss << "min " << min << ", max " << max;
	return ss.str();
}

}