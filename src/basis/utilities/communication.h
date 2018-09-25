
#ifndef SRC_BASIS_UTILITIES_COMMUNICATION_H_
#define SRC_BASIS_UTILITIES_COMMUNICATION_H_

#include "mpi.h"

#include <cstddef>
#include <vector>
#include <functional>

namespace espreso {

struct LoaderConfiguration;

class MPITools
{
	class Operations {
		friend class MPITools;
	public:
		struct TypedOperations {
			MPI_Op scan;
			MPI_Op max, min, sum;
		};

		MPI_Op mergeStatistics;
		TypedOperations SIZET, INT, LONG;

	private:
		Operations();
	public:
		Operations(Operations const&) = delete;
		void operator=(Operations const&) = delete;
	};

public:
	struct MPIType {
		MPI_Datatype type;
		size_t multiplier;

		MPIType(MPI_Datatype type): type(type), multiplier(1) {}
		MPIType(MPI_Datatype type, size_t multiplier): type(type), multiplier(multiplier) {}
	};

	template <typename Ttype>
	static MPITools::MPIType getType();

	static Operations& operations()
	{
		static Operations instance;
		return instance;
	}

	static Operations::TypedOperations& sizetOperations() { return operations().SIZET; }
	static Operations::TypedOperations& intOperations() { return operations().INT; }
	static Operations::TypedOperations& longOperations() { return operations().LONG; }
	static Operations::TypedOperations& eslocalOperations()
	{
			return sizeof(eslocal) == sizeof(int) ? operations().INT : operations().LONG;
	}

	struct MPIGroup {
		MPI_Comm communicator;
		int rank, size;

		MPIGroup();
	};

	class MPICommunicator {
		friend class MPITools;

	public:
		MPIGroup within, across;

		MPICommunicator();
	private:
		MPICommunicator(Operations const&) = delete;
		void operator=(MPICommunicator const&) = delete;
	};

	static MPICommunicator& nodes()
	{
		static MPICommunicator instance;
		return instance;
	}

	static MPIGroup& procs()
	{
		static MPIGroup instance;
		return instance;
	}


private:
	MPITools() = delete;
};

struct Communication {

	template <typename Ttype>
	static bool exchangeKnownSize(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbours, MPITools::MPIGroup &group = MPITools::procs());

	template <typename Ttype>
	static bool exchangeUnknownSize(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbours, MPITools::MPIGroup &group = MPITools::procs());

	template <typename Ttype>
	static bool exchangeUnknownSize(const std::vector<Ttype> &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbours, MPITools::MPIGroup &group = MPITools::procs());

	template <typename Ttype>
	static bool receiveLowerKnownSize(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbours, MPITools::MPIGroup &group = MPITools::procs());

	template <typename Ttype>
	static bool receiveUpperKnownSize(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbours, MPITools::MPIGroup &group = MPITools::procs());

	template <typename Ttype>
	static bool receiveUpperUnknownSize(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbours, MPITools::MPIGroup &group = MPITools::procs());

	template <typename Ttype>
	static bool gatherUnknownSize(const std::vector<Ttype> &sBuffer, std::vector<Ttype> &rBuffer, MPITools::MPIGroup &group = MPITools::procs());

	template <typename Ttype>
	static bool gatherUnknownSize(const std::vector<Ttype> &sBuffer, std::vector<Ttype> &rBuffer, std::vector<size_t> &offsets, MPITools::MPIGroup &group = MPITools::procs());

	template <typename Ttype>
	static bool allGatherUnknownSize(std::vector<Ttype> &data, MPITools::MPIGroup &group = MPITools::procs());

	template <typename Ttype>
	static bool broadcastUnknownSize(std::vector<Ttype> &buffer, MPITools::MPIGroup &group = MPITools::procs());

	template <typename Ttype>
	static bool balance(std::vector<Ttype> &buffer, const std::vector<size_t> &currentDistribution, const std::vector<size_t> &targetDistribution, MPITools::MPIGroup &group = MPITools::procs());

	template <typename Ttype>
	static bool allToAllV(const std::vector<Ttype> &sBuffer, std::vector<Ttype> &rBuffer, const std::vector<int> &ssize, const std::vector<int> &rsize, MPITools::MPIGroup &group = MPITools::procs());

	template <typename Ttype>
	static Ttype exscan(Ttype &value, MPITools::MPIGroup &group = MPITools::procs());

	template <typename Ttype>
	static std::vector<Ttype> getDistribution(Ttype size, MPITools::MPIGroup &group = MPITools::procs());

	template <typename Ttype>
	static bool sendVariousTargets(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &targets, MPITools::MPIGroup &group = MPITools::procs())
	{
		std::vector<int> sources;
		return sendVariousTargets(sBuffer, rBuffer, targets, sources, group);
	}

	template <typename Ttype>
	static bool sendVariousTargets(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &targets, std::vector<int> &sources, MPITools::MPIGroup &group = MPITools::procs());

	template <typename Ttype>
	static bool allToAllWithDataSizeAndTarget(const std::vector<Ttype> &sBuffer, std::vector<Ttype> &rBuffer, MPITools::MPIGroup &group = MPITools::procs());

	static void serialize(std::function<void(void)> fnc, MPITools::MPIGroup &group = MPITools::procs());

	static void createCommunicator(const LoaderConfiguration &configuration, MPITools::MPICommunicator &communicator);

private:
	template <typename Ttype>
	static Ttype exscan(Ttype &value, MPI_Op &operation, MPITools::MPIGroup &group);

	template <typename Ttype>
	static std::vector<Ttype> getDistribution(Ttype size, MPI_Op &operation, MPITools::MPIGroup &group);
};


}

#include "communication.hpp"




#endif /* SRC_BASIS_UTILITIES_COMMUNICATION_H_ */
