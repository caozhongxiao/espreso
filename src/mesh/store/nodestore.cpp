
#include "store.h"
#include "nodestore.h"

#include "../../basis/containers/point.h"
#include "../../basis/containers/serializededata.h"
#include "../../basis/utilities/utils.h"
#include "../../basis/utilities/communication.h"

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

  coordinates(NULL),
  ranks(NULL),

  idomains(NULL),
  ineighborOffsets(NULL),
  iranks(NULL),

  _dataBuffersPrepared(false)
{

}

size_t NodeStore::packedSize() const
{
	if (IDs == NULL || coordinates == NULL || idomains == NULL) {
		ESINFO(ERROR) << "ESPRESO internal error: invalid request for packedSize.";
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
			Esutils::packedSize(pintervals);
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
}

size_t NodeStore::packedDataSize() const
{
	size_t size = sizeof(size_t);
	for (size_t i = 0; i < data.size(); i++) {
		if (data[i]->names.size()) {
			size += data[i]->packedSize();
		}
	}
	return size;
}

void NodeStore::packData(char* &p) const
{
	size_t size = 0;
	for (size_t i = 0; i < data.size(); i++) {
		if (data[i]->names.size()) {
			++size;
		}
	}
	Esutils::pack(size, p);
	for (size_t i = 0; i < data.size(); i++) {
		if (data[i]->names.size()) {
			data[i]->pack(p);
		}
	}
}

void NodeStore::unpackData(const char* &p)
{
	size_t size;
	Esutils::unpack(size, p);
	for (size_t i = 0; i < size; i++) {
		if (i >= data.size()) {
			data.push_back(new NodeData());
		}
		data.back()->unpack(p);
	}
}

NodeStore::~NodeStore()
{
	if (IDs == NULL) { delete IDs; }
	if (elements == NULL) { delete elements; }

	if (coordinates == NULL) { delete coordinates; }
	if (ranks == NULL) { delete ranks; }

	if (idomains == NULL) { delete idomains; }
	if (ineighborOffsets == NULL) { delete ineighborOffsets; }
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
	Store::storedata(os, "ineighborOffsets", ineighborOffsets);
	Store::storedata(os, "iranks", iranks);
}

void NodeStore::permute(const std::vector<eslocal> &permutation, const std::vector<size_t> &distribution)
{
	this->distribution = distribution;

	if (IDs != NULL) { IDs->permute(permutation, distribution); }
	if (elements != NULL) { elements->permute(permutation, distribution); }

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

NodeData* NodeStore::appendData(const std::vector<std::string> &names)
{
	data.push_back(new NodeData(names));
	data.back()->decomposedData->resize(dintervals.size());
	for (size_t d = 0; d < dintervals.size(); d++) {
		data.back()->decomposedData[d].resize(dintervals[d].back().DOFOffset + dintervals[d].back().end - dintervals[d].back().begin);
	}
	return data.back();
}

NodeData* NodeStore::appendData(const std::vector<std::string> &names, std::vector<std::vector<double> > *data)
{
	this->data.push_back(new NodeData(names, data));
	return this->data.back();
}

size_t NodeData::packedSize() const
{
	size_t size = 0;
	for (size_t i = 0; i < names.size(); i++) {
		size += sizeof(size_t) + names[i].size();
	}
	return sizeof(size_t) + size + Esutils::packedSize(*gathredData);
}

void NodeData::pack(char* &p) const
{
	size_t size = names.size();
	Esutils::pack(size, p);
	for (size_t i = 0; i < names.size(); i++) {
		Esutils::pack(names[i], p);
	}
	Esutils::pack(*gathredData, p);
}

void NodeData::unpack(const char* &p)
{
	size_t size;
	Esutils::unpack(size, p);
	names.resize(size);
	for (size_t i = 0; i < names.size(); i++) {
		Esutils::unpack(names[i], p);
	}
	Esutils::unpack(*gathredData, p);
}

NodeData::NodeData()
: _delete(true)
{
	decomposedData = new std::vector<std::vector<double> >();
	gathredData = new std::vector<double>();
}

NodeData::NodeData(const std::vector<std::string> &names, std::vector<std::vector<double> > *data)
: names(names), decomposedData(data), _delete(false)
{
	if (decomposedData == NULL) {
		decomposedData = new std::vector<std::vector<double> >();
	}
	gathredData = new std::vector<double>();
}

NodeData::NodeData(NodeData &&other)
: names(std::move(other.names)), decomposedData(other.decomposedData), gathredData(other.gathredData),
  _delete(std::move(other._delete))
{
	other.decomposedData = NULL;
	other.gathredData = NULL;
	other._delete = false;
}

NodeData::NodeData(const NodeData &other)
: names(other.names), decomposedData(other.decomposedData),
  _delete(other._delete)
{
	if (_delete) {
		decomposedData = new std::vector<std::vector<double> >(*other.decomposedData);
	}
	gathredData = new std::vector<double>(*other.gathredData);
}

NodeData& NodeData::operator=(const NodeData &other)
{
	if (this != &other) {
		names = other.names;
		decomposedData = other.decomposedData;
		_delete = other._delete;
		if (_delete) {
			decomposedData = new std::vector<std::vector<double> >(*other.decomposedData);
		}
		gathredData = new std::vector<double>(*other.gathredData);
	}
	return *this;
}

NodeData::~NodeData()
{
	if (_delete && decomposedData != NULL) {
		delete decomposedData;
	}
	if (gathredData != NULL) {
		delete gathredData;
	}
}


