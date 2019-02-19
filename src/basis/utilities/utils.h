
#ifndef BASIS_UTILITIES_UTILS_H_
#define BASIS_UTILITIES_UTILS_H_

#include <vector>
#include <string>

namespace espreso {

struct Esutils
{
	template<typename Ttype>
	static Ttype sizesToOffsets(std::vector<Ttype> &sizes);

	template<typename Ttype>
	static std::vector<Ttype> sizesToOffsets(std::vector<std::vector<Ttype> > &sizes)
	{
		return sizesToOffsets(sizes, std::vector<Ttype>(sizes.front().size()));
	}

	template<typename Ttype>
	static std::vector<Ttype> sizesToOffsets(std::vector<std::vector<Ttype> > &sizes, const std::vector<Ttype> &offsets);

	template<typename Ttype, typename Tpermutation>
	static void permute(std::vector<Ttype> &data, const std::vector<Tpermutation> &permutation, size_t elementsize = 1);

	template<typename Ttype>
	static void threadDistributionToFullDistribution(std::vector<std::vector<Ttype> > &distribution);

	template<typename Ttype>
	static void threadDistributionToFullDistribution(std::vector<Ttype> &data, const std::vector<size_t> &distribution);

	template<typename Ttype>
	static void removeDuplicity(std::vector<Ttype> &data, size_t begin = 0);

	template<typename Ttype>
	static void sortAndRemoveDuplicity(std::vector<Ttype> &data, size_t begin = 0);

	template<typename Ttype>
	static void sortAndRemoveDuplicity(std::vector<std::vector<Ttype> > &data);

	template<typename Ttype>
	static void mergeThreadedUniqueData(std::vector<std::vector<Ttype> > &data);

	template<typename Ttype>
	static void mergeThreadedUniqueData(std::vector<std::vector<std::vector<Ttype> > > &data);

	template<typename Ttype>
	static void inplaceMerge(std::vector<Ttype> &data, const std::vector<size_t> &distribution);

	template<typename Ttype>
	static void inplaceMerge(std::vector<std::vector<Ttype> > &data);

	template<typename Ttype>
	static void sortWithInplaceMerge(std::vector<Ttype> &data, const std::vector<size_t> &distribution);

	template<typename Ttype>
	static void sortWithInplaceMerge(std::vector<std::vector<Ttype> > &data);

	template<typename Ttype>
	static void sortWithUniqueMerge(std::vector<std::vector<Ttype> > &data);

	template<typename Ttype>
	static void mergeAppendedData(std::vector<Ttype> &data, const std::vector<size_t> &distribution);

	template<typename Ttype>
	static typename std::vector<Ttype>::const_iterator max_element(const std::vector<Ttype> &elements);

	static std::string createDirectory(const std::vector<std::string> &path);
};

}

#include "utils.hpp"


#endif /* BASIS_UTILITIES_UTILS_H_ */
