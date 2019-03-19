
#ifndef SRC_BASIS_UTILITIES_COMMUNICATION_H_
#define SRC_BASIS_UTILITIES_COMMUNICATION_H_

#include "mpi.h"

#include <cstddef>
#include <string>
#include <vector>
#include <functional>

namespace espreso {

struct ProcessesReduction;

class MPIOperations {
	friend class MPITools;
public:
	struct TypedOperations {
		MPI_Op scan;
		MPI_Op max, min, sum;
	};

	MPI_Op mergeStatistics;
	TypedOperations SIZET, INT, LONG;

private:
	MPIOperations();
	~MPIOperations();
public:
	MPIOperations(MPIOperations const&) = delete;
	void operator=(MPIOperations const&) = delete;
};

struct MPIType {
	MPI_Datatype type;
	size_t multiplier;

	MPIType(MPI_Datatype type): type(type), multiplier(1) {}
	MPIType(MPI_Datatype type, size_t multiplier): type(type), multiplier(multiplier) {}
};

struct MPIGroup {
	MPI_Comm communicator;
	int rank, size;

	MPIGroup();
	MPIGroup(MPI_Comm &comm);
	~MPIGroup();
};

class MPISubset {
	friend class MPITools;

public:
	MPIGroup within, across, &origin;

	MPISubset(const ProcessesReduction &reduction, MPIGroup &origin);
private:
	MPISubset(MPISubset const&) = delete;
	void operator=(MPISubset const&) = delete;

	static void fillNodeColor();

	static int nodeRank;
};

class MPITools
{

public:
	static MPIOperations *operations;
	static MPIGroup *procs;
	static MPIGroup *instances;
	static MPIGroup *global;

	template <typename Ttype>
	static MPIType getType();

	static void init();
	static void destroy();

	static MPIOperations::TypedOperations& sizetOperations() { return operations->SIZET; }
	static MPIOperations::TypedOperations& intOperations() { return operations->INT; }
	static MPIOperations::TypedOperations& longOperations() { return operations->LONG; }
	static MPIOperations::TypedOperations& esintOperations()
	{
			return sizeof(esint) == sizeof(int) ? operations->INT : operations->LONG;
	}

private:
	MPITools() = delete;
};

struct Communication {

	template <typename Ttype>
	static bool exchangeKnownSize(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbours, MPIGroup *group = MPITools::procs);

	template <typename Ttype>
	static bool exchangeUnknownSize(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbours, MPIGroup *group = MPITools::procs);

	template <typename Ttype>
	static bool exchangeUnknownSize(const std::vector<Ttype> &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbours, MPIGroup *group = MPITools::procs);

	template <typename Ttype>
	static bool receiveLowerKnownSize(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbours, MPIGroup *group = MPITools::procs);

	template <typename Ttype>
	static bool receiveLowerUnknownSize(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbours, MPIGroup *group = MPITools::procs);

	template <typename Ttype>
	static bool receiveUpperKnownSize(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbours, MPIGroup *group = MPITools::procs);

	template <typename Ttype>
	static bool receiveUpperUnknownSize(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbours, MPIGroup *group = MPITools::procs);

	template <typename Ttype>
	static bool gatherUnknownSize(const std::vector<Ttype> &sBuffer, std::vector<Ttype> &rBuffer, MPIGroup *group = MPITools::procs);

	template <typename Ttype>
	static bool gatherUnknownSize(const std::vector<Ttype> &sBuffer, std::vector<Ttype> &rBuffer, std::vector<size_t> &offsets, MPIGroup *group = MPITools::procs);

	template <typename Ttype>
	static bool allGatherUnknownSize(std::vector<Ttype> &data, MPIGroup *group = MPITools::procs);

	template <typename Ttype>
	static bool broadcastUnknownSize(std::vector<Ttype> &buffer, MPIGroup *group = MPITools::procs);

	template <typename Ttype, typename Tdistribution>
	static bool balance(std::vector<Ttype> &buffer, const std::vector<Tdistribution> &currentDistribution, const std::vector<Tdistribution> &targetDistribution, MPIGroup *group = MPITools::procs);

	template <typename Ttype>
	static bool allToAllV(const std::vector<Ttype> &sBuffer, std::vector<Ttype> &rBuffer, const std::vector<int> &ssize, const std::vector<int> &rsize, MPIGroup *group = MPITools::procs);

	template <typename Ttype>
	static Ttype exscan(Ttype &value, MPIGroup *group = MPITools::procs);

	template <typename Ttype>
	static std::vector<Ttype> getDistribution(Ttype size, MPIGroup *group = MPITools::procs);

	template <typename Ttype>
	static bool sendVariousTargets(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &targets, MPIGroup *group = MPITools::procs)
	{
		std::vector<int> sources;
		return sendVariousTargets(sBuffer, rBuffer, targets, sources, group);
	}

	template <typename Ttype>
	static bool sendVariousTargets(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &targets, std::vector<int> &sources, MPIGroup *group = MPITools::procs);

	template <typename Ttype>
	static bool allToAllWithDataSizeAndTarget(const std::vector<Ttype> &sBuffer, std::vector<Ttype> &rBuffer, MPIGroup *group = MPITools::procs);

	static void serialize(std::function<void(void)> fnc, MPIGroup *group = MPITools::procs);

private:
	template <typename Ttype>
	static Ttype exscan(Ttype &value, MPI_Op &operation, MPIGroup *group);

	template <typename Ttype>
	static std::vector<Ttype> getDistribution(Ttype size, MPI_Op &operation, MPIGroup *group);
};


}

#include "communication.hpp"




#endif /* SRC_BASIS_UTILITIES_COMMUNICATION_H_ */
