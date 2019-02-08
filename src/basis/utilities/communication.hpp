
#include "communication.h"

#include <algorithm>
#include <cmath>

namespace espreso {

template <>
inline MPIType MPITools::getType<int>()
{
	return { MPI_INT };
}

template <>
inline MPIType MPITools::getType<uint>()
{
	return { MPI_UNSIGNED };
}

template <>
inline MPIType MPITools::getType<long>()
{
	return { MPI_LONG };
}

template <>
inline MPIType MPITools::getType<size_t>()
{
	return { MPI_UNSIGNED_LONG };
}

template <>
inline MPIType MPITools::getType<float>()
{
	return { MPI_FLOAT };
}

template <>
inline MPIType MPITools::getType<double>()
{
	return { MPI_DOUBLE, };
}

template <typename Ttype>
inline MPIType MPITools::getType()
{
	return { MPI_BYTE, sizeof(Ttype) };
}

template <typename Ttype>
bool Communication::exchangeKnownSize(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbours, MPIGroup *group)
{
	MPIType type = MPITools::getType<Ttype>();

	for (size_t n = 0; n < neighbours.size(); n++) {
		if (type.multiplier * sBuffer[n].size() > 1 << 30) {
			return false;
		}
		if (type.multiplier * rBuffer[n].size() > 1 << 30) {
			return false;
		}
	}
	std::vector<MPI_Request> req(2 * neighbours.size());

	for (size_t n = 0; n < neighbours.size(); n++) {
		// bullxmpi violate MPI standard (cast away constness)
		MPI_Isend(const_cast<Ttype*>(sBuffer[n].data()), type.multiplier * sBuffer[n].size(), type.type, neighbours[n], 0, group->communicator, req.data() + 2 * n);
	}

	for (size_t n = 0; n < neighbours.size(); n++) {
		// bullxmpi violate MPI standard (cast away constness)
		MPI_Irecv(const_cast<Ttype*>(rBuffer[n].data()), type.multiplier * rBuffer[n].size(), type.type, neighbours[n], 0, group->communicator, req.data() + 2 * n + 1);
	}

	MPI_Waitall(2 * neighbours.size(), req.data(), MPI_STATUSES_IGNORE);
	return true;
}


template <typename Ttype>
bool Communication::exchangeUnknownSize(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbours, MPIGroup *group)
{
	auto n2i = [ & ] (size_t neighbour) {
		return std::lower_bound(neighbours.begin(), neighbours.end(), neighbour) - neighbours.begin();
	};

	MPIType type = MPITools::getType<Ttype>();

	for (size_t n = 0; n < neighbours.size(); n++) {
		if (type.multiplier * sBuffer[n].size() > 1 << 30) {
			return false;
		}
	}

	std::vector<MPI_Request> req(neighbours.size());
	for (size_t n = 0; n < neighbours.size(); n++) {
		// bullxmpi violate MPI standard (cast away constness)
		MPI_Isend(const_cast<Ttype*>(sBuffer[n].data()), type.multiplier * sBuffer[n].size(), type.type, neighbours[n], 0, group->communicator, req.data() + n);
	}

	int flag;
	size_t counter = 0;
	MPI_Status status;
	while (counter < neighbours.size()) {
		MPI_Iprobe(MPI_ANY_SOURCE, 0, group->communicator, &flag, &status);
		if (flag) {
			int count;
			MPI_Get_count(&status, type.type, &count);
			rBuffer[n2i(status.MPI_SOURCE)].resize(count / type.multiplier);
			MPI_Recv(rBuffer[n2i(status.MPI_SOURCE)].data(), count, type.type, status.MPI_SOURCE, 0, group->communicator, MPI_STATUS_IGNORE);
			counter++;
		}
	}

	MPI_Waitall(neighbours.size(), req.data(), MPI_STATUSES_IGNORE);
	MPI_Barrier(group->communicator); // MPI_Iprobe(ANY_SOURCE) can be problem when calling this function more times
	return true;
}

template <typename Ttype>
bool Communication::exchangeUnknownSize(const std::vector<Ttype> &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbours, MPIGroup *group)
{
	auto n2i = [ & ] (size_t neighbour) {
		return std::lower_bound(neighbours.begin(), neighbours.end(), neighbour) - neighbours.begin();
	};

	MPIType type = MPITools::getType<Ttype>();
	if (type.multiplier * sBuffer.size() > 1 << 30) {
		return false;
	}

	std::vector<MPI_Request> req(neighbours.size());
	for (size_t n = 0; n < neighbours.size(); n++) {
		// bullxmpi violate MPI standard (cast away constness)
		MPI_Isend(const_cast<Ttype*>(sBuffer.data()), type.multiplier * sBuffer.size(), type.type, neighbours[n], 0, group->communicator, req.data() + n);
	}

	int flag;
	size_t counter = 0;
	MPI_Status status;
	while (counter < neighbours.size()) {
		MPI_Iprobe(MPI_ANY_SOURCE, 0, group->communicator, &flag, &status);
		if (flag) {
			int count;
			MPI_Get_count(&status, type.type, &count);
			rBuffer[n2i(status.MPI_SOURCE)].resize(count / type.multiplier);
			MPI_Recv(rBuffer[n2i(status.MPI_SOURCE)].data(), count, type.type, status.MPI_SOURCE, 0, group->communicator, MPI_STATUS_IGNORE);
			counter++;
		}
	}

	MPI_Waitall(neighbours.size(), req.data(), MPI_STATUSES_IGNORE);
	MPI_Barrier(group->communicator); // MPI_Iprobe(ANY_SOURCE) can be problem when calling this function more times
	return true;
}

template <typename Ttype>
bool Communication::receiveLowerKnownSize(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbours, MPIGroup *group)
{
	MPIType type = MPITools::getType<Ttype>();

	for (size_t n = 0; n < neighbours.size(); n++) {
		if (neighbours[n] > group->rank) {
			if (type.multiplier * sBuffer[n].size() > 1 << 30) {
				return false;
			}
		} else {
			if (type.multiplier * rBuffer[n].size() > 1 << 30) {
				return false;
			}
		}
	}

	std::vector<MPI_Request> req(neighbours.size());
	for (size_t n = 0; n < neighbours.size(); n++) {
		if (neighbours[n] > group->rank) {
			// bullxmpi violate MPI standard (cast away constness)
			MPI_Isend(const_cast<Ttype*>(sBuffer[n].data()), type.multiplier * sBuffer[n].size(), type.type, neighbours[n], 1, group->communicator, req.data() + n);
		}
		if (neighbours[n] < group->rank) {
			MPI_Irecv(rBuffer[n].data(), type.multiplier * rBuffer[n].size(), type.type, neighbours[n], 1, group->communicator, req.data() + n);
		}
	}

	MPI_Waitall(neighbours.size(), req.data(), MPI_STATUSES_IGNORE);
	return true;
}

template <typename Ttype>
bool Communication::receiveLowerUnknownSize(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbours, MPIGroup *group)
{
	MPIType type = MPITools::getType<Ttype>();
	for (size_t n = 0; n < neighbours.size(); n++) {
		if (type.multiplier * sBuffer[n].size() > 1 << 30) {
			return false;
		}
	}

	auto n2i = [ & ] (size_t neighbour) {
		return std::lower_bound(neighbours.begin(), neighbours.end(), neighbour) - neighbours.begin();
	};

	size_t rSize = 0;
	std::vector<MPI_Request> req(neighbours.size());
	for (size_t n = 0; n < neighbours.size(); n++) {
		if (group->rank < neighbours[n]) {
			// bullxmpi violate MPI standard (cast away constness)
			MPI_Isend(const_cast<Ttype*>(sBuffer[n].data()), type.multiplier * sBuffer[n].size(), type.type, neighbours[n], 0, group->communicator, req.data() + rSize++);
		}
	}

	int flag;
	size_t counter = neighbours.end() - std::lower_bound(neighbours.begin(), neighbours.end(), group->rank);
	MPI_Status status;
	while (counter < neighbours.size()) {
		MPI_Iprobe(MPI_ANY_SOURCE, 0, group->communicator, &flag, &status);
		if (flag) {
			int count;
			MPI_Get_count(&status, type.type, &count);
			rBuffer[n2i(status.MPI_SOURCE)].resize(count / type.multiplier);
			MPI_Recv(rBuffer[n2i(status.MPI_SOURCE)].data(), count, type.type, status.MPI_SOURCE, 0, group->communicator, MPI_STATUS_IGNORE);
			counter++;
		}
	}

	MPI_Waitall(rSize, req.data(), MPI_STATUSES_IGNORE);
	MPI_Barrier(group->communicator); // MPI_Iprobe(ANY_SOURCE) can be problem when calling this function more times
	return true;
}

template <typename Ttype>
bool Communication::receiveUpperKnownSize(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbours, MPIGroup *group)
{
	MPIType type = MPITools::getType<Ttype>();
	for (size_t n = 0; n < neighbours.size(); n++) {
		if (neighbours[n] < group->rank) {
			if (type.multiplier * sBuffer[n].size() > 1 << 30) {
				return false;
			}
		} else {
			if (type.multiplier * rBuffer[n].size() > 1 << 30) {
				return false;
			}
		}
	}

	std::vector<MPI_Request> req(neighbours.size());
	for (size_t n = 0; n < neighbours.size(); n++) {
		if (neighbours[n] < group->rank) {
			// bullxmpi violate MPI standard (cast away constness)
			MPI_Isend(const_cast<Ttype*>(sBuffer[n].data()), type.multiplier * sBuffer[n].size(), type.type, neighbours[n], 1, group->communicator, req.data() + n);
		}
		if (neighbours[n] > group->rank) {
			MPI_Irecv(rBuffer[n].data(), type.multiplier * rBuffer[n].size(), type.type, neighbours[n], 1, group->communicator, req.data() + n);
		}
	}

	MPI_Waitall(neighbours.size(), req.data(), MPI_STATUSES_IGNORE);
	return true;
}

template <typename Ttype>
bool Communication::receiveUpperUnknownSize(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbours, MPIGroup *group)
{
	MPIType type = MPITools::getType<Ttype>();
	for (size_t n = 0; n < neighbours.size() && neighbours[n] < group->rank; n++) {
		if (type.multiplier * sBuffer[n].size() > 1 << 30) {
			return false;
		}
	}

	auto n2i = [ & ] (size_t neighbour) {
		return std::lower_bound(neighbours.begin(), neighbours.end(), neighbour) - neighbours.begin();
	};

	size_t rSize = 0;
	std::vector<MPI_Request> req(neighbours.size());
	for (size_t n = 0; n < neighbours.size() && neighbours[n] < group->rank; n++) {
		// bullxmpi violate MPI standard (cast away constness)
		MPI_Isend(const_cast<Ttype*>(sBuffer[n].data()), type.multiplier * sBuffer[n].size(), type.type, neighbours[n], 0, group->communicator, req.data() + rSize++);
	}

	int flag;
	size_t counter = std::lower_bound(neighbours.begin(), neighbours.end(), group->rank) - neighbours.begin();
	MPI_Status status;
	while (counter < neighbours.size()) {
		MPI_Iprobe(MPI_ANY_SOURCE, 0, group->communicator, &flag, &status);
		if (flag) {
			int count;
			MPI_Get_count(&status, type.type, &count);
			rBuffer[n2i(status.MPI_SOURCE)].resize(count / type.multiplier);
			MPI_Recv(rBuffer[n2i(status.MPI_SOURCE)].data(), count, type.type, status.MPI_SOURCE, 0, group->communicator, MPI_STATUS_IGNORE);
			counter++;
		}
	}

	MPI_Waitall(rSize, req.data(), MPI_STATUSES_IGNORE);
	MPI_Barrier(group->communicator); // MPI_Iprobe(ANY_SOURCE) can be problem when calling this function more times
	return true;
}

template <typename Ttype>
bool Communication::gatherUnknownSize(const std::vector<Ttype> &sBuffer, std::vector<Ttype> &rBuffer, MPIGroup *group)
{
	std::vector<size_t> offsets;
	return gatherUnknownSize(sBuffer, rBuffer, offsets, group);
}

template <typename Ttype>
bool Communication::gatherUnknownSize(const std::vector<Ttype> &sBuffer, std::vector<Ttype> &rBuffer, std::vector<size_t> &offsets, MPIGroup *group)
{
	MPIType type = MPITools::getType<Ttype>();

	int size = type.multiplier * sBuffer.size();
	std::vector<int> rSizes(group->size), rOffsets(group->size);
	MPI_Gather(&size, 1, MPI_INT, rSizes.data(), 1, MPI_INT, 0, group->communicator);

	if (!group->rank) {
		size = 0;
		for (size_t i = 0; i < rSizes.size(); i++) {
			rOffsets[i] = size;
			size += rSizes[i];
		}
		rBuffer.resize(size / type.multiplier);
	}

	// bullxmpi violate MPI standard (cast away constness)
	MPI_Gatherv(const_cast<Ttype*>(sBuffer.data()), type.multiplier * sBuffer.size(), type.type, rBuffer.data(), rSizes.data(), rOffsets.data(), type.type, 0, group->communicator);

	offsets.resize(group->size);
	for (size_t i = 0; i < rOffsets.size(); i++) {
		offsets[i] = rOffsets[i] / type.multiplier;
	}
	return true;
}
template <typename Ttype>
bool Communication::allGatherUnknownSize(std::vector<Ttype> &data, MPIGroup *group)
{
	MPIType type = MPITools::getType<Ttype>();
	if (type.multiplier * data.size() > 1 << 30) {
		return false;
	}
	int size = type.multiplier * data.size();
	std::vector<int> rSizes(group->size), rOffsets(group->size);
	MPI_Allgather(&size, 1, MPI_INT, rSizes.data(), 1, MPI_INT, group->communicator);

	std::vector<Ttype> rdata;
	size = 0;
	for (size_t i = 0; i < rSizes.size(); i++) {
		rOffsets[i] = size;
		size += rSizes[i];
	}
	rdata.resize(size / type.multiplier);

	MPI_Allgatherv(data.data(), type.multiplier * data.size(), type.type, rdata.data(), rSizes.data(), rOffsets.data(), type.type, group->communicator);

	rdata.swap(data);
	return true;
}

template <typename Ttype>
bool Communication::broadcastUnknownSize(std::vector<Ttype> &buffer, MPIGroup *group)
{
	MPIType type = MPITools::getType<Ttype>();
	if (type.multiplier * buffer.size() > 1 << 30) {
		return false;
	}
	int size = buffer.size();
	MPI_Bcast(&size, 1, MPI_INT, 0, group->communicator);
	buffer.resize(size);
	MPI_Bcast(buffer.data(), type.multiplier * size, type.type, 0, group->communicator);
	return true;
}

template <typename Ttype, typename Tdistribution>
bool Communication::balance(std::vector<Ttype> &buffer, const std::vector<Tdistribution> &currentDistribution, const std::vector<Tdistribution> &targetDistribution, MPIGroup *group)
{
	MPIType type = MPITools::getType<Ttype>();
	if (type.multiplier * buffer.size() > 1 << 30) {
		return false;
	}

	std::vector<Ttype> result(targetDistribution[group->rank + 1] - targetDistribution[group->rank]);
	std::vector<int> ssize(group->size), sdisp(group->size), rsize(group->size), rdisp(group->size);

	auto fill = [&] (
			const std::vector<Tdistribution> &from, const std::vector<Tdistribution> &to,
			std::vector<int> &size, std::vector<int> &disp) {

		Tdistribution offset = 0;
		Tdistribution restSize = from[group->rank + 1] - from[group->rank];
		Tdistribution tIndex = std::lower_bound(to.begin(), to.end(), from[group->rank] + 1) - to.begin() - 1;
		Tdistribution tOffset = from[group->rank] - to[tIndex];
		while (restSize) {
			if (restSize < to[tIndex + 1] - to[tIndex] - tOffset) {
				size[tIndex] = restSize;
				disp[tIndex] = offset;
				restSize = 0;
			} else {
				size[tIndex] = to[tIndex + 1] - to[tIndex] - tOffset;
				disp[tIndex] = offset;
				restSize -= size[tIndex];
				offset += size[tIndex];
				++tIndex;
			}
			tOffset = 0;
		}

//		for (int r = 0; r < group->size; ++r) {
//			size[r] *= sizeof(Ttype);
//			disp[r] *= sizeof(Ttype);
//		}
	};

	fill(currentDistribution, targetDistribution, ssize, sdisp);
	fill(targetDistribution, currentDistribution, rsize, rdisp);

	std::vector<MPI_Request> requests(group->size + 1);
	int nrequests = 0;

	for (int r = 0; r < group->size; ++r) {
		if (rsize[r]) {
			MPI_Irecv(result.data() + rdisp[r], type.multiplier * rsize[r], type.type, r, 0, group->communicator, requests.data() + nrequests++);
		}
	}

	for (int r = 0; r < group->size; ++r) {
		if (ssize[r]) {
			MPI_Isend(buffer.data() + sdisp[r], type.multiplier * ssize[r], type.type, r, 0, group->communicator, requests.data() + nrequests++);
		}
	}

	MPI_Waitall(nrequests, requests.data(), MPI_STATUSES_IGNORE);
//	MPI_Alltoallv(buffer.data(), ssize.data(), sdisp.data(), MPI_BYTE, result.data(), rsize.data(), rdisp.data(), MPI_BYTE, group->communicator);
	buffer.swap(result);

	return true;
}

template <typename Ttype>
bool Communication::allToAllV(const std::vector<Ttype> &sBuffer, std::vector<Ttype> &rBuffer, const std::vector<int> &ssize, const std::vector<int> &rsize, MPIGroup *group)
{
	MPIType type = MPITools::getType<Ttype>();
	if (type.multiplier * sBuffer.size() > 1 << 30) {
		return false;
	}
	std::vector<int> _ssize = ssize, _rsize = rsize;
	std::vector<int> sdisp(group->size), rdisp(group->size);
	for (int r = 0; r < group->size; r++) {
		_ssize[r] *= type.multiplier;
		_rsize[r] *= type.multiplier;
	}
	for (int r = 1; r < group->size; r++) {
		sdisp[r] = sdisp[r - 1] + _ssize[r - 1];
		rdisp[r] = rdisp[r - 1] + _rsize[r - 1];
	}
	MPI_Alltoallv(sBuffer.data(), _ssize.data(), sdisp.data(), type.type, rBuffer.data(), _rsize.data(), rdisp.data(), type.type, group->communicator);
	return true;
}

template <>
inline size_t Communication::exscan(size_t &value, MPIGroup *group)
{
	return exscan(value, MPITools::sizetOperations().scan, group);
}

template <>
inline int Communication::exscan(int &value, MPIGroup *group)
{
	return exscan(value, MPITools::intOperations().scan, group);
}

template <>
inline long Communication::exscan(long &value, MPIGroup *group)
{
	return exscan(value, MPITools::longOperations().scan, group);
}

template <typename Ttype>
Ttype Communication::exscan(Ttype &value, MPI_Op &operation, MPIGroup *group)
{
	Ttype size = value;
	if (group->size == 1) {
		value = 0;
		return size;
	}

	MPI_Exscan(&size, &value, sizeof(Ttype), MPI_BYTE, operation, group->communicator);

	size = value + size;
	MPI_Bcast(&size, sizeof(Ttype), MPI_BYTE, group->size - 1, group->communicator);
	if (group->rank == 0) {
		value = 0;
	}
	MPI_Barrier(group->communicator);

	return size;
}

template <>
inline std::vector<size_t> Communication::getDistribution(size_t size, MPIGroup *group)
{
	return getDistribution(size, MPITools::sizetOperations().scan, group);
}

template <>
inline std::vector<int> Communication::getDistribution(int size, MPIGroup *group)
{
	return getDistribution(size, MPITools::intOperations().scan, group);
}

template <>
inline std::vector<long> Communication::getDistribution(long size, MPIGroup *group)
{
	return getDistribution(size, MPITools::longOperations().scan, group);
}

template <typename Ttype>
std::vector<Ttype> Communication::getDistribution(Ttype size, MPI_Op &operation, MPIGroup *group)
{
	std::vector<Ttype> result(group->size + 1);
	Ttype esize = size;
	Communication::exscan(esize, operation, group);

	MPI_Allgather(&esize, sizeof(Ttype), MPI_BYTE, result.data(), sizeof(Ttype), MPI_BYTE, group->communicator);
	result.back() = esize + size;
	MPI_Bcast(&result.back(), sizeof(Ttype), MPI_BYTE, group->size - 1, group->communicator);
	return result;
}

template <typename Ttype>
bool Communication::sendVariousTargets(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &targets, std::vector<int> &sources, MPIGroup *group)
{
	MPIType type = MPITools::getType<Ttype>();
	for (size_t n = 0; n < targets.size(); n++) {
		if (type.multiplier * sBuffer[n].size() > 1 << 30) {
			return false;
		}
	}

	std::vector<int> smsgcounter(group->size);
	std::vector<int> rmsgcounter(group->size);
	for (size_t n = 0; n < targets.size(); n++) {
		smsgcounter[targets[n]] = 1;
	}

	MPI_Allreduce(smsgcounter.data(), rmsgcounter.data(), group->size, MPI_INT, MPI_SUM, group->communicator);

	std::vector<MPI_Request> req(targets.size());
	for (size_t t = 0; t < targets.size(); t++) {
		MPI_Isend(const_cast<Ttype*>(sBuffer[t].data()), type.multiplier * sBuffer[t].size(), type.type, targets[t], 0, group->communicator, req.data() + t);
	}

	int flag;
	int counter = 0;
	MPI_Status status;
	sources.clear();
	std::vector<std::vector<Ttype> > tmpBuffer;
	tmpBuffer.reserve(rmsgcounter[group->rank]);
	while (counter < rmsgcounter[group->rank]) {
		MPI_Iprobe(MPI_ANY_SOURCE, 0, group->communicator, &flag, &status);
		if (flag) {
			int count;
			MPI_Get_count(&status, type.type, &count);
			tmpBuffer.push_back(std::vector<Ttype>(count / type.multiplier));
			MPI_Recv(tmpBuffer.back().data(), count, type.type, status.MPI_SOURCE, 0, group->communicator, MPI_STATUS_IGNORE);
			sources.push_back(status.MPI_SOURCE);
			counter++;
		}
	}

	std::vector<int> permutation(sources.size());
	for (size_t i = 0; i < sources.size(); i++) {
		permutation[i] = i;
	}
	std::sort(permutation.begin(), permutation.end(), [&] (int i, int j) { return sources[i] < sources[j]; });
	rBuffer.resize(tmpBuffer.size());
	for (size_t i = 0; i < permutation.size(); i++) {
		rBuffer[i].swap(tmpBuffer[permutation[i]]);
	}

	std::sort(sources.begin(), sources.end());

	MPI_Waitall(targets.size(), req.data(), MPI_STATUSES_IGNORE);
	MPI_Barrier(group->communicator); // MPI_Iprobe(ANY_SOURCE) can be problem when calling this function more times
	return true;
}

template <typename Ttype>
bool Communication::allToAllWithDataSizeAndTarget(const std::vector<Ttype> &sBuffer, std::vector<Ttype> &rBuffer, MPIGroup *group)
{
	MPIType type = MPITools::getType<Ttype>();
	size_t levels = std::ceil(std::log2(group->size));

	std::vector<Ttype> prevsend, send, recv;
	recv.reserve(sBuffer.size());

	send = sBuffer;

	MPI_Status status;
	int recvsize, recvmidsize;

	auto movebefore = [&] (std::vector<Ttype> &data, int rank, size_t begin, size_t end) {
		size_t pos = begin;
		while (pos < end && (int)data[pos + 1] < rank) {
			pos += (size_t)data[pos];
		}
		return pos;
	};

	size_t mybegin = movebefore(send, group->rank, 0, send.size());
	size_t myend = movebefore(send, group->rank + 1, mybegin, send.size());
	rBuffer.insert(rBuffer.end(), send.begin() + mybegin, send.begin() + myend);

	int left = 0, right = group->size, mid;
	for (size_t l = 0; l < levels && left + 1 < right; l++) {
		mid = left + (right - left) / 2 + (right - left) % 2;
		if (group->rank < mid) {
			// LOWER half to UPPER half
			if (group->rank + (mid - left) >= right) {
				// PRE :
				// send: l1, l2, l3, ME, u1, u2, u3

				// POST:
				// SEND: ME -> u1
				// RECV:

				// send: l1, l2, l3
				size_t my = movebefore(send, group->rank, 0, send.size());
				size_t upper = movebefore(send, mid, my, send.size());

				if (type.multiplier * (send.size() - upper) > 1 << 30) {
					return false;
				}
				MPI_Send(send.data() + upper, type.multiplier * (send.size() - upper), type.type, mid, 0, group->communicator);
				send.resize(my);
			} else {
				// PRE :
				// send: l1, l2(ME), l3, l4, u1, u2, u3


				// POST:
				// SEND: u1, u2, u3 -> u2;
				// RECV: l1, l2, l3, l4 <- u2

				// sBuffer: l1, l2, l3, l4, u1, u2, u3
				// send: l1, l1, l2, l2, l3, l3, l4, l4

				size_t upper = movebefore(send, mid, 0, send.size());
				if (type.multiplier * (send.size() - upper) > 1 << 30) {
					return false;
				}
				MPI_Send(send.data() + upper, type.multiplier * (send.size() - upper), type.type, group->rank + (mid - left), 0, group->communicator);

				MPI_Probe(group->rank + (mid - left), 0, group->communicator, &status);
				MPI_Get_count(&status, type.type, &recvsize);
				recv.resize(recvsize / type.multiplier);
				MPI_Recv(recv.data(), recvsize, type.type, group->rank + (mid - left), 0, group->communicator, MPI_STATUS_IGNORE);

				send.swap(prevsend);
				send.clear();

				size_t recvbegin = 0;
				size_t recvend = recvbegin;
				size_t sendbegin = 0;
				size_t sendend = sendbegin;
				for (int r = left; r < mid; r++) {
					recvbegin = recvend;
					recvend = movebefore(recv, r + 1, recvbegin, recv.size());
					sendbegin = sendend;
					sendend = movebefore(prevsend, r + 1, sendbegin, prevsend.size());
					if (r == group->rank) {
						rBuffer.insert(rBuffer.end(), recv.begin() + recvbegin, recv.begin() + recvend);
					} else {
						send.insert(send.end(), recv.begin() + recvbegin, recv.begin() + recvend);
						send.insert(send.end(), prevsend.begin() + sendbegin, prevsend.begin() + sendend);
					}
				}
			}
			right = mid;
		} else {
			// UPPER half to LOWER half

			size_t upper = movebefore(send, mid, 0, send.size());

			MPI_Probe(group->rank - (mid - left), 0, group->communicator, &status);
			MPI_Get_count(&status, type.type, &recvsize);
			recv.resize(recvsize / type.multiplier);
			MPI_Recv(recv.data(), recvsize, type.type, group->rank - (mid - left), 0, group->communicator, MPI_STATUS_IGNORE);
			if (type.multiplier * upper > 1 << 30) {
				return false;
			}
			MPI_Send(send.data(), type.multiplier * upper, type.type, group->rank - (mid - left), 0, group->communicator);

			recvmidsize = recvsize / type.multiplier;
			if (mid - left > right - mid && group->rank == mid) {
				// l1, l2, l3, l4, u1(ME), u2, u3
				// RECV: l4
				MPI_Probe(mid - 1, 0, group->communicator, &status);
				MPI_Get_count(&status, type.type, &recvsize);
				recv.resize(recv.size() + recvsize / type.multiplier);
				MPI_Recv(recv.data() + recvmidsize, recvsize, type.type, mid - 1, 0, group->communicator, MPI_STATUS_IGNORE);
			}

			// PRE :
			// send: l1, l2, l3, l4, u1(ME), u2, u3

			// POST:
			// SEND: l1, l2, l3, l4 -> l1;
			// RECV: u1, u2, u3 <- l1
			// RECV: u1, u2, u3 <- l4 (recvsize > recvmidsize)

			// sBuffer: l1, l2, l3, l4, u1, u2, u3
			// send: l1, l1, l2, l2, l3, l3, l4, l4

			send.swap(prevsend);
			send.clear();

			size_t recvbegin = 0;
			size_t recvend = recvbegin;
			size_t recvmidbegin = recvmidsize;
			size_t recvmidend = recvmidbegin;
			size_t sendbegin = movebefore(prevsend, mid, 0, prevsend.size());
			size_t sendend = sendbegin;
			for (int r = mid; r < right; r++) {
				recvbegin = recvend;
				recvend = movebefore(recv, r + 1, recvbegin, recvmidsize);
				recvmidbegin = recvmidend;
				recvmidend = movebefore(recv, r + 1, recvmidbegin, recv.size());
				sendbegin = sendend;
				sendend = movebefore(prevsend, r + 1, sendbegin, prevsend.size());
				if (r == group->rank) {
					rBuffer.insert(rBuffer.end(), recv.begin() + recvbegin, recv.begin() + recvend);
					rBuffer.insert(rBuffer.end(), recv.begin() + recvmidbegin, recv.begin() + recvmidend);
				} else {
					send.insert(send.end(), recv.begin() + recvbegin, recv.begin() + recvend);
					send.insert(send.end(), recv.begin() + recvmidbegin, recv.begin() + recvmidend);
					send.insert(send.end(), prevsend.begin() + sendbegin, prevsend.begin() + sendend);
				}
			}

			left = mid;
		}
	}
	return true;
}

inline void Communication::serialize(std::function<void(void)> fnc, MPIGroup *group)
{
	for (int r = 0; r < group->size; ++r) {
		if (r == group->rank) {
			fnc();
		}
		MPI_Barrier(group->communicator);
	}
	MPI_Barrier(group->communicator);
}

}



