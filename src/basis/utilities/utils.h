
#ifndef BASIS_UTILITIES_UTILS_H_
#define BASIS_UTILITIES_UTILS_H_

#include <vector>

namespace espreso {
namespace utils {

	template<typename Ttype>
	Ttype sizesToOffsets(std::vector<Ttype> &sizes);

	template<typename Ttype>
	std::vector<Ttype> sizesToOffsets(std::vector<std::vector<Ttype> > &sizes, const std::vector<Ttype> &offsets);

	template<typename Ttype>
	std::vector<Ttype> sizesToOffsets(std::vector<std::vector<Ttype> > &sizes)
	{
		return sizesToOffsets(sizes, std::vector<Ttype>(sizes.front().size()));
	}

	template<typename Ttype, typename Tpermutation>
	void permute(std::vector<Ttype> &data, const std::vector<Tpermutation> &permutation, size_t elementsize = 1);

	template<typename Ttype>
	void threadDistributionToFullDistribution(std::vector<std::vector<Ttype> > &distribution);

	template<typename Ttype>
	void threadDistributionToFullDistribution(std::vector<Ttype> &data, const std::vector<size_t> &distribution);

	template<typename Ttype>
	void removeDuplicity(std::vector<Ttype> &data, size_t begin = 0);

	template<typename Ttype>
	void sortAndRemoveDuplicity(std::vector<Ttype> &data, size_t begin = 0);

	template<typename Ttype>
	void sortAndRemoveDuplicity(std::vector<std::vector<Ttype> > &data);

	template<typename Ttype>
	void mergeThreadedUniqueData(std::vector<std::vector<Ttype> > &data);

	template<typename Ttype>
	void mergeThreadedUniqueData(std::vector<std::vector<std::vector<Ttype> > > &data);

	template<typename Ttype>
	void inplaceMerge(std::vector<Ttype> &data, const std::vector<size_t> &distribution);

	template<typename Ttype>
	void inplaceMerge(std::vector<std::vector<Ttype> > &data);

	template<typename Ttype>
	void sortWithInplaceMerge(std::vector<Ttype> &data, const std::vector<size_t> &distribution);

	template<typename Ttype>
	void sortWithInplaceMerge(std::vector<std::vector<Ttype> > &data);

	template<typename Ttype>
	void sortWithUniqueMerge(std::vector<std::vector<Ttype> > &data);

	template<typename Ttype>
	void mergeAppendedData(std::vector<Ttype> &data, const std::vector<size_t> &distribution);
};

}

#include "utils.hpp"


#endif /* BASIS_UTILITIES_UTILS_H_ */
