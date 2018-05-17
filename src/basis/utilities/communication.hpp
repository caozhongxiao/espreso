
#include "communication.h"

#include <algorithm>
#include <cmath>

namespace espreso {

template <typename Ttype>
bool Communication::exchangeKnownSize(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbours)
{
	std::vector<MPI_Request> req(2 * neighbours.size());
	for (size_t n = 0; n < neighbours.size(); n++) {
		// bullxmpi violate MPI standard (cast away constness)
		MPI_Isend(const_cast<Ttype*>(sBuffer[n].data()), sizeof(Ttype) * sBuffer[n].size(), MPI_BYTE, neighbours[n], 0, environment->MPICommunicator, req.data() + 2 * n);
	}

	for (size_t n = 0; n < neighbours.size(); n++) {
		// bullxmpi violate MPI standard (cast away constness)
		MPI_Irecv(const_cast<Ttype*>(rBuffer[n].data()), sizeof(Ttype) * rBuffer[n].size(), MPI_BYTE, neighbours[n], 0, environment->MPICommunicator, req.data() + 2 * n + 1);
	}

	MPI_Waitall(2 * neighbours.size(), req.data(), MPI_STATUSES_IGNORE);
	return true;
}


template <typename Ttype>
bool Communication::exchangeUnknownSize(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbours)
{
	auto n2i = [ & ] (size_t neighbour) {
		return std::lower_bound(neighbours.begin(), neighbours.end(), neighbour) - neighbours.begin();
	};

	std::vector<MPI_Request> req(neighbours.size());
	for (size_t n = 0; n < neighbours.size(); n++) {
		// bullxmpi violate MPI standard (cast away constness)
		MPI_Isend(const_cast<Ttype*>(sBuffer[n].data()), sizeof(Ttype) * sBuffer[n].size(), MPI_BYTE, neighbours[n], 0, environment->MPICommunicator, req.data() + n);
	}

	int flag;
	size_t counter = 0;
	MPI_Status status;
	while (counter < neighbours.size()) {
		MPI_Iprobe(MPI_ANY_SOURCE, 0, environment->MPICommunicator, &flag, &status);
		if (flag) {
			int count;
			MPI_Get_count(&status, MPI_BYTE, &count);
			rBuffer[n2i(status.MPI_SOURCE)].resize(count / sizeof(Ttype));
			MPI_Recv(rBuffer[n2i(status.MPI_SOURCE)].data(), count, MPI_BYTE, status.MPI_SOURCE, 0, environment->MPICommunicator, MPI_STATUS_IGNORE);
			counter++;
		}
	}

	MPI_Waitall(neighbours.size(), req.data(), MPI_STATUSES_IGNORE);
	MPI_Barrier(environment->MPICommunicator); // MPI_Iprobe(ANY_SOURCE) can be problem when calling this function more times
	return true;
}

template <typename Ttype>
bool Communication::exchangeUnknownSize(const std::vector<Ttype> &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbours)
{
	auto n2i = [ & ] (size_t neighbour) {
		return std::lower_bound(neighbours.begin(), neighbours.end(), neighbour) - neighbours.begin();
	};

	std::vector<MPI_Request> req(neighbours.size());
	for (size_t n = 0; n < neighbours.size(); n++) {
		// bullxmpi violate MPI standard (cast away constness)
		MPI_Isend(const_cast<Ttype*>(sBuffer.data()), sizeof(Ttype) * sBuffer.size(), MPI_BYTE, neighbours[n], 0, environment->MPICommunicator, req.data() + n);
	}

	int flag;
	size_t counter = 0;
	MPI_Status status;
	while (counter < neighbours.size()) {
		MPI_Iprobe(MPI_ANY_SOURCE, 0, environment->MPICommunicator, &flag, &status);
		if (flag) {
			int count;
			MPI_Get_count(&status, MPI_BYTE, &count);
			rBuffer[n2i(status.MPI_SOURCE)].resize(count / sizeof(Ttype));
			MPI_Recv(rBuffer[n2i(status.MPI_SOURCE)].data(), count, MPI_BYTE, status.MPI_SOURCE, 0, environment->MPICommunicator, MPI_STATUS_IGNORE);
			counter++;
		}
	}

	MPI_Waitall(neighbours.size(), req.data(), MPI_STATUSES_IGNORE);
	MPI_Barrier(environment->MPICommunicator); // MPI_Iprobe(ANY_SOURCE) can be problem when calling this function more times
	return true;
}

template <typename Ttype>
bool Communication::receiveLowerKnownSize(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbours)
{
	std::vector<MPI_Request> req(neighbours.size());
	for (size_t n = 0; n < neighbours.size(); n++) {
		if (neighbours[n] > environment->MPIrank) {
			// bullxmpi violate MPI standard (cast away constness)
			MPI_Isend(const_cast<Ttype*>(sBuffer[n].data()), sizeof(Ttype) * sBuffer[n].size(), MPI_BYTE, neighbours[n], 1, environment->MPICommunicator, req.data() + n);
		}
		if (neighbours[n] < environment->MPIrank) {
			MPI_Irecv(rBuffer[n].data(), sizeof(Ttype) * rBuffer[n].size(), MPI_BYTE, neighbours[n], 1, environment->MPICommunicator, req.data() + n);
		}
	}

	MPI_Waitall(neighbours.size(), req.data(), MPI_STATUSES_IGNORE);
	return true;
}

template <typename Ttype>
bool Communication::receiveUpperKnownSize(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbours)
{
	std::vector<MPI_Request> req(neighbours.size());
	for (size_t n = 0; n < neighbours.size(); n++) {
		if (neighbours[n] < environment->MPIrank) {
			// bullxmpi violate MPI standard (cast away constness)
			MPI_Isend(const_cast<Ttype*>(sBuffer[n].data()), sizeof(Ttype) * sBuffer[n].size(), MPI_BYTE, neighbours[n], 1, environment->MPICommunicator, req.data() + n);
		}
		if (neighbours[n] > environment->MPIrank) {
			MPI_Irecv(rBuffer[n].data(), sizeof(Ttype) * rBuffer[n].size(), MPI_BYTE, neighbours[n], 1, environment->MPICommunicator, req.data() + n);
		}
	}

	MPI_Waitall(neighbours.size(), req.data(), MPI_STATUSES_IGNORE);
	return true;
}

template <typename Ttype>
bool Communication::receiveUpperUnknownSize(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbours)
{
	auto n2i = [ & ] (size_t neighbour) {
		return std::lower_bound(neighbours.begin(), neighbours.end(), neighbour) - neighbours.begin();
	};

	size_t rSize = 0;
	std::vector<MPI_Request> req(neighbours.size());
	for (size_t n = 0; n < neighbours.size() && neighbours[n] < environment->MPIrank; n++) {
		// bullxmpi violate MPI standard (cast away constness)
		MPI_Isend(const_cast<Ttype*>(sBuffer[n].data()), sizeof(Ttype) * sBuffer[n].size(), MPI_BYTE, neighbours[n], 0, environment->MPICommunicator, req.data() + rSize++);
	}

	int flag;
	size_t counter = std::lower_bound(neighbours.begin(), neighbours.end(), environment->MPIrank) - neighbours.begin();
	MPI_Status status;
	while (counter < neighbours.size()) {
		MPI_Iprobe(MPI_ANY_SOURCE, 0, environment->MPICommunicator, &flag, &status);
		if (flag) {
			int count;
			MPI_Get_count(&status, MPI_BYTE, &count);
			rBuffer[n2i(status.MPI_SOURCE)].resize(count / sizeof(Ttype));
			MPI_Recv(rBuffer[n2i(status.MPI_SOURCE)].data(), count, MPI_BYTE, status.MPI_SOURCE, 0, environment->MPICommunicator, MPI_STATUS_IGNORE);
			counter++;
		}
	}

	MPI_Waitall(rSize, req.data(), MPI_STATUSES_IGNORE);
	MPI_Barrier(environment->MPICommunicator); // MPI_Iprobe(ANY_SOURCE) can be problem when calling this function more times
	return true;
}

template <typename Ttype>
bool Communication::gatherUnknownSize(const std::vector<Ttype> &sBuffer, std::vector<Ttype> &rBuffer)
{
	std::vector<size_t> offsets;
	return gatherUnknownSize(sBuffer, rBuffer, offsets);
}

template <typename Ttype>
bool Communication::gatherUnknownSize(const std::vector<Ttype> &sBuffer, std::vector<Ttype> &rBuffer, std::vector<size_t> &offsets)
{
	int size = sizeof(Ttype) * sBuffer.size();
	std::vector<int> rSizes(environment->MPIsize), rOffsets(environment->MPIsize);
	MPI_Gather(&size, sizeof(int), MPI_BYTE, rSizes.data(), sizeof(int), MPI_BYTE, 0, environment->MPICommunicator);

	if (!environment->MPIrank) {
		size = 0;
		for (size_t i = 0; i < rSizes.size(); i++) {
			rOffsets[i] = size;
			size += rSizes[i];
		}
		rBuffer.resize(size / sizeof(Ttype));
	}

	// bullxmpi violate MPI standard (cast away constness)
	MPI_Gatherv(const_cast<Ttype*>(sBuffer.data()), sBuffer.size() * sizeof(Ttype), MPI_BYTE, rBuffer.data(), rSizes.data(), rOffsets.data(), MPI_BYTE, 0, environment->MPICommunicator);

	offsets.resize(environment->MPIsize);
	for (size_t i = 0; i < rOffsets.size(); i++) {
		offsets[i] = rOffsets[i] / sizeof(Ttype);
	}
	return true;
}
template <typename Ttype>
bool Communication::allGatherUnknownSize(std::vector<Ttype> &data)
{
	int size = sizeof(Ttype) * data.size();
	std::vector<int> rSizes(environment->MPIsize), rOffsets(environment->MPIsize);
	MPI_Allgather(&size, sizeof(int), MPI_BYTE, rSizes.data(), sizeof(int), MPI_BYTE, environment->MPICommunicator);

	std::vector<Ttype> rdata;
	size = 0;
	for (size_t i = 0; i < rSizes.size(); i++) {
		rOffsets[i] = size;
		size += rSizes[i];
	}
	rdata.resize(size / sizeof(Ttype));

	MPI_Allgatherv(data.data(), data.size() * sizeof(Ttype), MPI_BYTE, rdata.data(), rSizes.data(), rOffsets.data(), MPI_BYTE, environment->MPICommunicator);

	rdata.swap(data);
	return true;
}

template <typename Ttype>
bool Communication::broadcastUnknownSize(std::vector<Ttype> &buffer)
{
	int size = buffer.size();
	MPI_Bcast(&size, sizeof(int), MPI_BYTE, 0, environment->MPICommunicator);
	buffer.resize(size);
	MPI_Bcast(buffer.data(), sizeof(Ttype) * size, MPI_BYTE, 0, environment->MPICommunicator);
	return true;
}

template <typename Ttype>
bool Communication::balance(std::vector<Ttype> &buffer, const std::vector<size_t> &currentDistribution, const std::vector<size_t> &targetDistribution)
{
	std::vector<Ttype> result(targetDistribution[environment->MPIrank + 1] - targetDistribution[environment->MPIrank]);
	std::vector<int> ssize(environment->MPIsize), sdisp(environment->MPIsize), rsize(environment->MPIsize), rdisp(environment->MPIsize);

	auto fill = [] (
			const std::vector<size_t> &from, const std::vector<size_t> &to,
			std::vector<int> &size, std::vector<int> &disp) {

		size_t offset = 0;
		size_t restSize = from[environment->MPIrank + 1] - from[environment->MPIrank];
		size_t tIndex = std::lower_bound(to.begin(), to.end(), from[environment->MPIrank] + 1) - to.begin() - 1;
		size_t tOffset = from[environment->MPIrank] - to[tIndex];
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

//		for (int r = 0; r < environment->MPIsize; ++r) {
//			size[r] *= sizeof(Ttype);
//			disp[r] *= sizeof(Ttype);
//		}
	};

	fill(currentDistribution, targetDistribution, ssize, sdisp);
	fill(targetDistribution, currentDistribution, rsize, rdisp);

	std::vector<MPI_Request> requests(environment->MPIsize + 1);
	int nrequests = 0;

	for (int r = 0; r < environment->MPIsize; ++r) {
		if (rsize[r]) {
			MPI_Irecv(result.data() + rdisp[r], rsize[r] * sizeof(Ttype), MPI_BYTE, r, 0, environment->MPICommunicator, requests.data() + nrequests++);
		}
	}

	for (int r = 0; r < environment->MPIsize; ++r) {
		if (ssize[r]) {
			MPI_Isend(buffer.data() + sdisp[r], ssize[r] * sizeof(Ttype), MPI_BYTE, r, 0, environment->MPICommunicator, requests.data() + nrequests++);
		}
	}

	MPI_Waitall(nrequests, requests.data(), MPI_STATUSES_IGNORE);
//	MPI_Alltoallv(buffer.data(), ssize.data(), sdisp.data(), MPI_BYTE, result.data(), rsize.data(), rdisp.data(), MPI_BYTE, environment->MPICommunicator);
	buffer.swap(result);

	return true;
}

template <typename Ttype>
bool Communication::allToAllV(const std::vector<Ttype> &sBuffer, std::vector<Ttype> &rBuffer, const std::vector<int> &ssize, const std::vector<int> &rsize)
{
	std::vector<int> _ssize = ssize, _rsize = rsize;
	std::vector<int> sdisp(environment->MPIsize), rdisp(environment->MPIsize);
	for (int r = 0; r < environment->MPIsize; r++) {
		_ssize[r] *= sizeof(Ttype);
		_rsize[r] *= sizeof(Ttype);
	}
	for (int r = 1; r < environment->MPIsize; r++) {
		sdisp[r] = sdisp[r - 1] + _ssize[r - 1];
		rdisp[r] = rdisp[r - 1] + _rsize[r - 1];
	}
	MPI_Alltoallv(sBuffer.data(), _ssize.data(), sdisp.data(), MPI_BYTE, rBuffer.data(), _rsize.data(), rdisp.data(), MPI_BYTE, environment->MPICommunicator);
	return true;
}

template <typename Ttype>
Ttype Communication::exscan(Ttype &value, MPI_Op &operation)
{
	Ttype size = value;
	if (environment->MPIsize == 1) {
		value = 0;
		return size;
	}

	MPI_Exscan(&size, &value, sizeof(Ttype), MPI_BYTE, operation, environment->MPICommunicator);

	size = value + size;
	MPI_Bcast(&size, sizeof(Ttype), MPI_BYTE, environment->MPIsize - 1, environment->MPICommunicator);
	if (environment->MPIrank == 0) {
		value = 0;
	}
	MPI_Barrier(environment->MPICommunicator);

	return size;
}

template <typename Ttype>
std::vector<Ttype> Communication::getDistribution(Ttype size, MPI_Op operation)
{
	std::vector<Ttype> result(environment->MPIsize + 1);
	Ttype esize = size;
	Communication::exscan(esize, operation);

	MPI_Allgather(&esize, sizeof(Ttype), MPI_BYTE, result.data(), sizeof(Ttype), MPI_BYTE, environment->MPICommunicator);
	result.back() = esize + size;
	MPI_Bcast(&result.back(), sizeof(Ttype), MPI_BYTE, environment->MPIsize - 1, environment->MPICommunicator);
	return result;
}

template <typename Ttype>
bool Communication::sendVariousTargets(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &targets, std::vector<int> &sources)
{
	std::vector<int> smsgcounter(environment->MPIsize);
	std::vector<int> rmsgcounter(environment->MPIsize);
	for (size_t n = 0; n < targets.size(); n++) {
		smsgcounter[targets[n]] = 1;
	}

	MPI_Allreduce(smsgcounter.data(), rmsgcounter.data(), environment->MPIsize, MPI_INT, MPI_SUM, environment->MPICommunicator);

	std::vector<MPI_Request> req(targets.size());
	for (size_t t = 0; t < targets.size(); t++) {
		MPI_Isend(const_cast<Ttype*>(sBuffer[t].data()), sizeof(Ttype) * sBuffer[t].size(), MPI_BYTE, targets[t], 0, environment->MPICommunicator, req.data() + t);
	}

	int flag;
	int counter = 0;
	MPI_Status status;
	sources.clear();
	std::vector<std::vector<Ttype> > tmpBuffer;
	tmpBuffer.reserve(rmsgcounter[environment->MPIrank]);
	while (counter < rmsgcounter[environment->MPIrank]) {
		MPI_Iprobe(MPI_ANY_SOURCE, 0, environment->MPICommunicator, &flag, &status);
		if (flag) {
			int count;
			MPI_Get_count(&status, MPI_BYTE, &count);
			tmpBuffer.push_back(std::vector<Ttype>(count / sizeof(Ttype)));
			MPI_Recv(tmpBuffer.back().data(), count, MPI_BYTE, status.MPI_SOURCE, 0, environment->MPICommunicator, MPI_STATUS_IGNORE);
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
	MPI_Barrier(environment->MPICommunicator); // MPI_Iprobe(ANY_SOURCE) can be problem when calling this function more times
	return true;
}

template <typename Ttype>
bool Communication::allToAllWithDataSizeAndTarget(const std::vector<Ttype> &sBuffer, std::vector<Ttype> &rBuffer)
{
	size_t levels = std::ceil(std::log2(environment->MPIsize));

	std::vector<Ttype> prevsend, send, recv;
	recv.reserve(sBuffer.size());

	send = sBuffer;

	MPI_Status status;
	int recvsize, recvmidsize;

	auto movebefore = [&] (std::vector<Ttype> &data, int rank, size_t begin, size_t end) {
		size_t pos = begin;
		while (pos < end && data[pos + 1] < rank) {
			pos += data[pos];
		}
		return pos;
	};

	size_t mybegin = movebefore(send, environment->MPIrank, 0, send.size());
	size_t myend = movebefore(send, environment->MPIrank + 1, mybegin, send.size());
	rBuffer.insert(rBuffer.end(), send.begin() + mybegin, send.begin() + myend);

	int left = 0, right = environment->MPIsize, mid;
	for (size_t l = 0; l < levels && left + 1 < right; l++) {
		mid = left + (right - left) / 2 + (right - left) % 2;
		if (environment->MPIrank < mid) {
			// LOWER half to UPPER half
			if (environment->MPIrank + (mid - left) >= right) {
				// PRE :
				// send: l1, l2, l3, ME, u1, u2, u3

				// POST:
				// SEND: ME -> u1
				// RECV:

				// send: l1, l2, l3
				size_t my = movebefore(send, environment->MPIrank, 0, send.size());
				size_t upper = movebefore(send, mid, my, send.size());

				MPI_Send(send.data() + upper, send.size() - upper, MPI_INT, mid, 0, environment->MPICommunicator);
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
				MPI_Send(send.data() + upper, send.size() - upper, MPI_INT, environment->MPIrank + (mid - left), 0, environment->MPICommunicator);

				MPI_Probe(environment->MPIrank + (mid - left), 0, environment->MPICommunicator, &status);
				MPI_Get_count(&status, MPI_INT, &recvsize);
				recv.resize(recvsize);
				MPI_Recv(recv.data(), recvsize, MPI_INT, environment->MPIrank + (mid - left), 0, environment->MPICommunicator, MPI_STATUS_IGNORE);

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
					if (r == environment->MPIrank) {
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

			MPI_Probe(environment->MPIrank - (mid - left), 0, environment->MPICommunicator, &status);
			MPI_Get_count(&status, MPI_INT, &recvsize);
			recv.resize(recvsize);
			MPI_Recv(recv.data(), recvsize, MPI_INT, environment->MPIrank - (mid - left), 0, environment->MPICommunicator, MPI_STATUS_IGNORE);
			MPI_Send(send.data(), upper, MPI_INT, environment->MPIrank - (mid - left), 0, environment->MPICommunicator);

			recvmidsize = recvsize;
			if (mid - left > right - mid && environment->MPIrank == mid) {
				// l1, l2, l3, l4, u1(ME), u2, u3
				// RECV: l4
				MPI_Probe(mid - 1, 0, environment->MPICommunicator, &status);
				MPI_Get_count(&status, MPI_INT, &recvsize);
				recv.resize(recv.size() + recvsize);
				MPI_Recv(recv.data() + recvmidsize, recvsize, MPI_INT, mid - 1, 0, environment->MPICommunicator, MPI_STATUS_IGNORE);
				recvsize += recvmidsize;
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
				if (r == environment->MPIrank) {
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

inline void Communication::serialize(std::function<void(void)> fnc)
{
	for (int r = 0; r < environment->MPIsize; ++r) {
		if (r == environment->MPIrank) {
			fnc();
		}
		MPI_Barrier(environment->MPICommunicator);
	}
	MPI_Barrier(environment->MPICommunicator);
}

}



