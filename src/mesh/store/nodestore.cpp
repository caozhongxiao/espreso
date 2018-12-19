
#include "store.h"
#include "nodestore.h"
#include "statisticsstore.h"

#include "../mesh.h"

#include "../../basis/containers/point.h"
#include "../../basis/containers/serializededata.h"
#include "../../basis/utilities/utils.h"
#include "../../basis/utilities/communication.h"

#include "../../globals/run.h"

#include "../../config/ecf/environment.h"

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

size_t NodeStore::packedSize() const
{
	if (IDs == NULL || coordinates == NULL || idomains == NULL) {
		ESINFO(ERROR) << "ESPRESO internal error: invalid request for packedSize.";
	}

	size_t datasize = sizeof(size_t);
	for (size_t i = 0; i < data.size(); i++) {
		if (data[i]->names.size()) {
			datasize += sizeof(int) + Esutils::packedSize(data[i]->names);
		}
	}

	return
			Esutils::packedSize(size) +
			Esutils::packedSize(uniqueOffset) +
			Esutils::packedSize(uniqueSize) +
			Esutils::packedSize(uniqueTotalSize) +
			IDs->packedSize() +
			coordinates->packedSize() +
			idomains->packedSize() +
			Esutils::packedSize(externalIntervals) +
			Esutils::packedSize(pintervals) +
			Esutils::packedSize(dintervals) +
			datasize +
			Esutils::packedSize(dcenter) +
			Esutils::packedSize(center);
}

void NodeStore::pack(char* &p) const
{
	Esutils::pack(size, p);
	Esutils::pack(uniqueOffset, p);
	Esutils::pack(uniqueSize, p);
	Esutils::pack(uniqueTotalSize, p);
	IDs->pack(p);
	coordinates->pack(p);
	idomains->pack(p);
	Esutils::pack(externalIntervals, p);
	Esutils::pack(pintervals, p);
	Esutils::pack(dintervals, p);
	Esutils::pack(dcenter, p);
	Esutils::pack(center, p);

	size_t size = 0;
	for (size_t i = 0; i < data.size(); i++) {
		if (data[i]->names.size()) {
			size += 1;
		}
	}
	Esutils::pack(size, p);
	for (size_t i = 0; i < data.size(); i++) {
		if (data[i]->names.size()) {
			Esutils::pack(data[i]->dimension, p);
			Esutils::pack(data[i]->names, p);
		}
	}
}

void NodeStore::unpack(const char* &p)
{
	if (IDs == NULL) {
		IDs = new serializededata<eslocal, eslocal>(1, tarray<eslocal>(1, 0));
	}
	if (coordinates == NULL) {
		coordinates = new serializededata<eslocal, Point>(1, tarray<Point>(1, 0));
	}
	if (idomains == NULL) {
		idomains = new serializededata<eslocal, eslocal>(tarray<eslocal>(1, 0), tarray<eslocal>(1, 0));
	}

	Esutils::unpack(size, p);
	Esutils::unpack(uniqueOffset, p);
	Esutils::unpack(uniqueSize, p);
	Esutils::unpack(uniqueTotalSize, p);
	IDs->unpack(p);
	coordinates->unpack(p);
	idomains->unpack(p);
	Esutils::unpack(externalIntervals, p);
	Esutils::unpack(pintervals, p);
	Esutils::unpack(dintervals, p);
	Esutils::unpack(dcenter, p);
	Esutils::unpack(center, p);

	int dimension;
	size_t datasize;
	Esutils::unpack(datasize, p);
	for (size_t i = 0; i < datasize; i++) {
		if (i >= data.size()) {
			Esutils::unpack(dimension, p);
			data.push_back(new NodeData(dimension, {}));
		}
		Esutils::unpack(data[i]->names, p);
	}
}

size_t NodeStore::packedDataSize() const
{
	size_t size = 0;
	for (size_t i = 0; i < data.size(); i++) {
		if (data[i]->names.size()) {
			size += Esutils::packedSize(data[i]->data);
		}
	}
	return size;
}

void NodeStore::packData(char* &p) const
{
	for (size_t i = 0; i < data.size(); i++) {
		if (data[i]->names.size()) {
			Esutils::pack(data[i]->data, p);
		}
	}
}

void NodeStore::unpackData(const char* &p)
{
	for (size_t i = 0; i < data.size(); i++) {
		Esutils::unpack(data[i]->data, p);
	}
}

NodeStore::~NodeStore()
{
	if (IDs == NULL) { delete IDs; }
	if (elements == NULL) { delete elements; }

	if (originCoordinates == NULL) { delete originCoordinates; }
	if (coordinates == NULL) { delete coordinates; }
	if (ranks == NULL) { delete ranks; }

	if (idomains == NULL) { delete idomains; }
	if (iranks == NULL) { delete iranks; }
}

void NodeStore::store(const std::string &file)
{
	std::ofstream os(file + std::to_string(environment->MPIrank) + ".txt");

	Store::storedata(os, "IDs", IDs);
	Store::storedata(os, "elements", elements);

	Store::storedata(os, "coordinates", coordinates);
	Store::storedata(os, "ranks", ranks);

	Store::storedata(os, "idomains", idomains);
	Store::storedata(os, "iranks", iranks);
}

void NodeStore::permute(const std::vector<eslocal> &permutation, const std::vector<size_t> &distribution)
{
	this->distribution = distribution;

	if (IDs != NULL) { IDs->permute(permutation, distribution); }
	if (elements != NULL) { elements->permute(permutation, distribution); }

	if (originCoordinates != NULL) { originCoordinates->permute(permutation, distribution); }
	if (coordinates != NULL) { coordinates->permute(permutation, distribution); }
	if (ranks != NULL) { ranks->permute(permutation, distribution); }
}

std::vector<eslocal> NodeStore::gatherNodeDistribution()
{
	return Store::gatherDistribution(size);
}

std::vector<eslocal> NodeStore::gatherUniqueNodeDistribution()
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

void NodeData::statistics(const tarray<eslocal> &nodes, eslocal totalsize, Statistics *statistics)
{
	for (int d = 0; d <= names.size(); d++) {
		(statistics + d)->reset();
	}

	auto nranks = run::mesh->nodes->ranks->begin();
	if (names.size() == 1) {
		eslocal prev = 0;
		for (auto n = nodes.begin(); n != nodes.end(); prev = *n++) {
			nranks += *n - prev;
			if (*nranks->begin() == environment->MPIrank) {
				statistics->min = std::min(statistics->min, data[*n * dimension]);
				statistics->max = std::max(statistics->max, data[*n * dimension]);
				statistics->avg += data[*n * dimension];
				statistics->norm += data[*n * dimension] * data[*n * dimension];
			}
		}
	} else {
		eslocal prev = 0;
		for (auto n = nodes.begin(); n != nodes.end(); prev = *n++) {
			nranks += *n - prev;
			if (*nranks->begin() == environment->MPIrank) {
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
	MPI_Allreduce(statistics, global.data(), sizeof(Statistics) * names.size(), MPI_BYTE, MPITools::operations().mergeStatistics, environment->MPICommunicator);
	memcpy(statistics, global.data(), sizeof(Statistics) * names.size());

	for (size_t i = 0; i < names.size(); i++) {
		(statistics + i)->avg /= totalsize;
		(statistics + i)->norm = std::sqrt((statistics + i)->norm);
	}
}


