
#include "store.h"
#include "nodestore.h"
#include "statisticsstore.h"

#include "mesh/mesh.h"

#include "esinfo/mpiinfo.h"
#include "esinfo/meshinfo.h"
#include "basis/containers/point.h"
#include "basis/containers/serializededata.h"
#include "basis/utilities/packing.h"

using namespace espreso;

NodeStore::NodeStore()
: size(0),
  uniqueOffset(0),
  uniqueSize(0),
  uniqueTotalSize(0),

  distribution({0, 0}),

  IDs(NULL),
  elements(NULL),

  originCoordinates(NULL),
  coordinates(NULL),
  ranks(NULL),

  idomains(NULL),
  iranks(NULL)
{

}

size_t NodeStore::packedFullSize() const
{
	size_t packedSize = 0;

	packedSize += utils::packedSize(size);
	packedSize += utils::packedSize(uniqueOffset);
	packedSize += utils::packedSize(uniqueSize);
	packedSize += utils::packedSize(uniqueTotalSize);
	packedSize += utils::packedSize(distribution);

	packedSize += utils::packedSize(IDs);
	packedSize += utils::packedSize(elements);
	packedSize += utils::packedSize(originCoordinates);
	packedSize += utils::packedSize(coordinates);
	packedSize += utils::packedSize(ranks);

	packedSize += utils::packedSize(externalIntervals);
	packedSize += utils::packedSize(pintervals);
	packedSize += utils::packedSize(dintervals);

	packedSize += utils::packedSize(idomains);
	packedSize += utils::packedSize(iranks);

	packedSize += utils::packedSize(center);
	packedSize += utils::packedSize(dcenter);
	packedSize += utils::packedSize(min);
	packedSize += utils::packedSize(max);
	packedSize += utils::packedSize(lmin);
	packedSize += utils::packedSize(lmax);

	packedSize += utils::packedSize(data.size());
	for (size_t i = 0; i < data.size(); i++) {
		packedSize += data[i]->packedSize();
	}

	return packedSize;
}

void NodeStore::packFull(char* &p) const
{
	utils::pack(size, p);
	utils::pack(uniqueOffset, p);
	utils::pack(uniqueSize, p);
	utils::pack(uniqueTotalSize, p);
	utils::pack(distribution, p);

	utils::pack(IDs, p);
	utils::pack(elements, p);
	utils::pack(originCoordinates, p);
	utils::pack(coordinates, p);
	utils::pack(ranks, p);

	utils::pack(externalIntervals, p);
	utils::pack(pintervals, p);
	utils::pack(dintervals, p);

	utils::pack(idomains, p);
	utils::pack(iranks, p);

	utils::pack(center, p);
	utils::pack(dcenter, p);
	utils::pack(min, p);
	utils::pack(max, p);
	utils::pack(lmin, p);
	utils::pack(lmax, p);

	utils::pack(data.size(), p);
	for (size_t i = 0; i < data.size(); i++) {
		data[i]->pack(p);
	}
}

void NodeStore::unpackFull(const char* &p)
{
	utils::unpack(size, p);
	utils::unpack(uniqueOffset, p);
	utils::unpack(uniqueSize, p);
	utils::unpack(uniqueTotalSize, p);
	utils::unpack(distribution, p);

	utils::unpack(IDs, p);
	utils::unpack(elements, p);
	utils::unpack(originCoordinates, p);
	utils::unpack(coordinates, p);
	utils::unpack(ranks, p);

	utils::unpack(externalIntervals, p);
	utils::unpack(pintervals, p);
	utils::unpack(dintervals, p);

	utils::unpack(idomains, p);
	utils::unpack(iranks, p);

	utils::unpack(center, p);
	utils::unpack(dcenter, p);
	utils::unpack(min, p);
	utils::unpack(max, p);
	utils::unpack(lmin, p);
	utils::unpack(lmax, p);

	size_t size;
	utils::unpack(size, p);
	for (size_t i = 0; i < size; i++) {
		data.push_back(new NodeData(p));
	}
}

size_t NodeStore::packedSize() const
{
	size_t datasize = sizeof(size_t);
	for (size_t i = 0; i < data.size(); i++) {
		if (data[i]->names.size()) {
			datasize += sizeof(int) + utils::packedSize(data[i]->names);
		}
	}

	return
			utils::packedSize(size) +
			utils::packedSize(uniqueOffset) +
			utils::packedSize(uniqueSize) +
			utils::packedSize(uniqueTotalSize) +
			IDs->packedSize() +
			coordinates->packedSize() +
			idomains->packedSize() +
			utils::packedSize(externalIntervals) +
			utils::packedSize(pintervals) +
			utils::packedSize(dintervals) +
			datasize +
			utils::packedSize(dcenter) +
			utils::packedSize(center);
}

void NodeStore::pack(char* &p) const
{
	utils::pack(size, p);
	utils::pack(uniqueOffset, p);
	utils::pack(uniqueSize, p);
	utils::pack(uniqueTotalSize, p);
	IDs->pack(p);
	coordinates->pack(p);
	idomains->pack(p);
	utils::pack(externalIntervals, p);
	utils::pack(pintervals, p);
	utils::pack(dintervals, p);
	utils::pack(dcenter, p);
	utils::pack(center, p);

	size_t size = 0;
	for (size_t i = 0; i < data.size(); i++) {
		if (data[i]->names.size()) {
			size += 1;
		}
	}
	utils::pack(size, p);
	for (size_t i = 0; i < data.size(); i++) {
		if (data[i]->names.size()) {
			utils::pack(data[i]->dimension, p);
			utils::pack(data[i]->names, p);
		}
	}
}

void NodeStore::unpack(const char* &p)
{
	if (IDs == NULL) {
		IDs = new serializededata<esint, esint>(1, tarray<esint>(1, 0));
	}
	if (coordinates == NULL) {
		coordinates = new serializededata<esint, Point>(1, tarray<Point>(1, 0));
	}
	if (idomains == NULL) {
		idomains = new serializededata<esint, esint>(tarray<esint>(1, 0), tarray<esint>(1, 0));
	}

	utils::unpack(size, p);
	utils::unpack(uniqueOffset, p);
	utils::unpack(uniqueSize, p);
	utils::unpack(uniqueTotalSize, p);
	IDs->unpack(p);
	coordinates->unpack(p);
	idomains->unpack(p);
	utils::unpack(externalIntervals, p);
	utils::unpack(pintervals, p);
	utils::unpack(dintervals, p);
	utils::unpack(dcenter, p);
	utils::unpack(center, p);

	int dimension;
	size_t datasize;
	utils::unpack(datasize, p);
	for (size_t i = 0; i < datasize; i++) {
		if (i >= data.size()) {
			utils::unpack(dimension, p);
			data.push_back(new NodeData(dimension, {}));
		}
		utils::unpack(data[i]->names, p);
	}
}

size_t NodeStore::packedDataSize() const
{
	size_t size = 0;
	for (size_t i = 0; i < data.size(); i++) {
		if (data[i]->names.size()) {
			size += utils::packedSize(data[i]->data);
		}
	}
	return size;
}

void NodeStore::packData(char* &p) const
{
	for (size_t i = 0; i < data.size(); i++) {
		if (data[i]->names.size()) {
			utils::pack(data[i]->data, p);
		}
	}
}

void NodeStore::unpackData(const char* &p)
{
	for (size_t i = 0; i < data.size(); i++) {
		utils::unpack(data[i]->data, p);
	}
}

NodeStore::~NodeStore()
{
	if (IDs != NULL) { delete IDs; }
	if (elements != NULL) { delete elements; }

	if (originCoordinates != NULL) { delete originCoordinates; }
	if (coordinates != NULL) { delete coordinates; }
	if (ranks != NULL) { delete ranks; }

	if (idomains != NULL) { delete idomains; }
	if (iranks != NULL) { delete iranks; }

	for (size_t i = 0; i < data.size(); i++) {
		delete data[i];
	}
}

void NodeStore::store(const std::string &file)
{
	std::ofstream os(file + std::to_string(info::mpi::rank) + ".txt");

	Store::storedata(os, "IDs", IDs);
	Store::storedata(os, "elements", elements);

	Store::storedata(os, "coordinates", coordinates);
	Store::storedata(os, "ranks", ranks);

	Store::storedata(os, "idomains", idomains);
	Store::storedata(os, "iranks", iranks);
}

void NodeStore::permute(const std::vector<esint> &permutation, const std::vector<size_t> &distribution)
{
	this->distribution = distribution;

	if (IDs != NULL) { IDs->permute(permutation, distribution); }
	if (elements != NULL) { elements->permute(permutation, distribution); }

	if (originCoordinates != NULL) { originCoordinates->permute(permutation, distribution); }
	if (coordinates != NULL) { coordinates->permute(permutation, distribution); }
	if (ranks != NULL) { ranks->permute(permutation, distribution); }
}

std::vector<esint> NodeStore::gatherNodeDistribution()
{
	return Store::gatherDistribution(size);
}

std::vector<esint> NodeStore::gatherUniqueNodeDistribution()
{
	return Store::gatherDistribution(uniqueSize);
}

NodeData* NodeStore::appendData(int dimension, const std::vector<std::string> &names)
{
	data.push_back(new NodeData(dimension, names));
	data.back()->data.resize(dimension * size);
	return data.back();
}

NodeData::NodeData(int dimension, const std::vector<std::string> &names)
: dimension(dimension), names(names)
{

}

NodeData::NodeData(const char* &packedData)
{
	utils::unpack(dimension, packedData);
	utils::unpack(names, packedData);
	utils::unpack(data, packedData);
}

size_t NodeData::packedSize()
{
	size_t packetSize = 0;
	packetSize += utils::packedSize(dimension);
	packetSize += utils::packedSize(names);
	packetSize += utils::packedSize(data);
	return packetSize;
}

void NodeData::pack(char *&p)
{
	utils::pack(dimension, p);
	utils::pack(names, p);
	utils::pack(data, p);
}

void NodeData::statistics(const tarray<esint> &nodes, esint totalsize, Statistics *statistics) const
{
	for (size_t d = 0; d < names.size(); d++) {
		(statistics + d)->reset();
	}

	auto nranks = info::mesh->nodes->ranks->begin();
	if (names.size() == 1) {
		esint prev = 0;
		for (auto n = nodes.begin(); n != nodes.end(); prev = *n++) {
			nranks += *n - prev;
			if (*nranks->begin() == info::mpi::rank) {
				statistics->min = std::min(statistics->min, data[*n * dimension]);
				statistics->max = std::max(statistics->max, data[*n * dimension]);
				statistics->avg += data[*n * dimension];
				statistics->norm += data[*n * dimension] * data[*n * dimension];
			}
		}
	} else {
		esint prev = 0;
		for (auto n = nodes.begin(); n != nodes.end(); prev = *n++) {
			nranks += *n - prev;
			if (*nranks->begin() == info::mpi::rank) {
				double value = 0;
				for (int d = 0; d < dimension; d++) {
					value += data[*n * dimension + d] * data[*n * dimension + d];
					(statistics + d + 1)->min = std::min((statistics + d + 1)->min, data[*n * dimension + d]);
					(statistics + d + 1)->max = std::max((statistics + d + 1)->max, data[*n * dimension + d]);
					(statistics + d + 1)->avg += data[*n * dimension + d];
					(statistics + d + 1)->norm += data[*n * dimension + d] * data[*n * dimension + d];
				}
				value = std::sqrt(value);
				statistics->min = std::min(statistics->min, value);
				statistics->max = std::max(statistics->max, value);
				statistics->avg += value;
				statistics->norm += value * value;
			}
		}
	}

	std::vector<Statistics> global(names.size());
	MPI_Allreduce(statistics, global.data(), sizeof(Statistics) * names.size(), MPI_BYTE, MPITools::operations->mergeStatistics, info::mpi::comm);
	memcpy(statistics, global.data(), sizeof(Statistics) * names.size());

	for (size_t i = 0; i < names.size(); i++) {
		(statistics + i)->avg /= totalsize;
		(statistics + i)->norm = std::sqrt((statistics + i)->norm);
	}
}

double NodeData::norm() const
{
	esint foreignPrefix = info::mesh->nodes->size - info::mesh->nodes->uniqueSize;
	double square = 0;
	for (size_t i = foreignPrefix * dimension; i < data.size(); ++i) {
		square += data[i] * data[i];
	}

	double sum = 0;
	MPI_Allreduce(&square, &sum, 1, MPI_DOUBLE, MPI_SUM, info::mpi::comm);

	return std::sqrt(sum);
}

double NodeData::maxabs() const
{
	esint foreignPrefix = info::mesh->nodes->size - info::mesh->nodes->uniqueSize;
	double max = 0;
	for (size_t i = foreignPrefix * dimension; i < data.size(); ++i) {
		max = std::max(std::fabs(data[i]), max);
	}

	double gmax = 0;
	MPI_Allreduce(&max, &gmax, 1, MPI_DOUBLE, MPI_MAX, info::mpi::comm);

	return gmax;
}

