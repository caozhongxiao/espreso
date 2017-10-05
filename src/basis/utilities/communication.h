
#ifndef SRC_BASIS_UTILITIES_COMMUNICATION_H_
#define SRC_BASIS_UTILITIES_COMMUNICATION_H_

namespace espreso {

struct Communication {

	template <typename Ttype>
	static bool exchangeKnownSize(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbours);

	template <typename Ttype>
	static bool exchangeUnknownSize(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbours);

	template <typename Ttype>
	static bool receiveLowerKnownSize(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbours);

	template <typename Ttype>
	static bool receiveUpperKnownSize(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbours);

	template <typename Ttype>
	static bool receiveUpperUnknownSize(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbours);

	template <typename Ttype>
	static bool gatherUnknownSize(const std::vector<Ttype> &sBuffer, std::vector<Ttype> &rBuffer);

	template <typename Ttype>
	static bool gatherUnknownSize(const std::vector<Ttype> &sBuffer, std::vector<Ttype> &rBuffer, std::vector<size_t> &offsets);

	template <typename Ttype>
	static bool broadcastUnknownSize(std::vector<Ttype> &buffer);

	template <typename Ttype>
	static Ttype exscan(Ttype &value);

	template <typename Ttype>
	static bool sendVariousTargets(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &targets)
	{
		std::vector<int> sources;
		return sendVariousTargets(sBuffer, rBuffer, targets, sources);
	}

	template <typename Ttype>
	static bool sendVariousTargets(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &targets, std::vector<int> &sources);
};


}

#include "communication.hpp"




#endif /* SRC_BASIS_UTILITIES_COMMUNICATION_H_ */
