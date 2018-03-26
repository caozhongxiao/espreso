
#include "store.h"
#include "elementstore.h"

#include "../elements/element.h"

#include "../../basis/containers/serializededata.h"
#include "../../config/ecf/environment.h"

using namespace espreso;

ElementStore::ElementStore(std::vector<Element*> &eclasses)
: size(0),
  distribution({0, 0}),

  IDs(NULL),
  nodes(NULL),

  body(NULL),
  material(NULL),
  regions(NULL),
  epointers(NULL),

  neighbors(NULL),

  firstDomain(0),
  ndomains(1),
  nclusters(1),

  regionMaskSize(1),
  ecounters(static_cast<int>(Element::CODE::SIZE)),

  _eclasses(eclasses)
{

}

size_t ElementStore::packedSize() const
{
	if (nodes == NULL || epointers == NULL) {
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
			nodes->packedSize() +
			sizeof(size_t) + epointers->datatarray().size() * sizeof(int) +
			Esutils::packedSize(firstDomain) +
			Esutils::packedSize(ndomains) +
			Esutils::packedSize(nclusters) +
			Esutils::packedSize(clusters) +
			Esutils::packedSize(elementsDistribution) +
			Esutils::packedSize(ecounters) +
			Esutils::packedSize(eintervals) +
			datasize;
}

void ElementStore::pack(char* &p) const
{
	Esutils::pack(size, p);
	nodes->pack(p);
	if (epointers != NULL) {
		std::vector<int> eindices;
		eindices.reserve(epointers->datatarray().size());

		size_t threads = environment->OMP_NUM_THREADS;
		for (size_t t = 0; t < threads; t++) {
			for (size_t i = this->distribution[t]; i < this->distribution[t + 1]; ++i) {
				eindices.push_back(epointers->datatarray()[i] - _eclasses[t]);
			}
		}
		Esutils::pack(eindices, p);
	}
	Esutils::pack(firstDomain, p);
	Esutils::pack(ndomains, p);
	Esutils::pack(nclusters, p);
	Esutils::pack(clusters, p);
	Esutils::pack(elementsDistribution, p);
	Esutils::pack(ecounters, p);
	Esutils::pack(eintervals, p);

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

void ElementStore::unpack(const char* &p)
{
	if (nodes == NULL) {
		nodes = new serializededata<eslocal, eslocal>(tarray<eslocal>(1, 0), tarray<eslocal>(1, 0));
	}
	if (epointers == NULL) {
		epointers = new serializededata<eslocal, Element*>(1, tarray<Element*>(1, 0));
	}

	Esutils::unpack(size, p);
	nodes->unpack(p);
	if (epointers != NULL) {
		std::vector<int> eindices;
		Esutils::unpack(eindices, p);
		if (epointers != NULL) {
			delete epointers;
		}
		epointers = new serializededata<eslocal, Element*>(1, tarray<Element*>(1, size));
		for (size_t i = 0; i < size; ++i) {
			epointers->datatarray()[i] = &_eclasses[0][eindices[i]];
		}
	}
	Esutils::unpack(firstDomain, p);
	Esutils::unpack(ndomains, p);
	Esutils::unpack(nclusters, p);
	Esutils::unpack(clusters, p);
	Esutils::unpack(elementsDistribution, p);
	Esutils::unpack(ecounters, p);
	Esutils::unpack(eintervals, p);

	int dimension;
	size_t datasize;
	Esutils::unpack(datasize, p);
	for (size_t i = 0; i < datasize; i++) {
		if (i >= data.size()) {
			Esutils::unpack(dimension, p);
			data.push_back(new ElementData(dimension));
		}
		Esutils::unpack(data[i]->names, p);
	}
}

size_t ElementStore::packedDataSize() const
{
	size_t size = 0;
	for (size_t i = 0; i < data.size(); i++) {
		if (data[i]->names.size()) {
			size += Esutils::packedSize(*data[i]->data);
		}
	}
	return size;
}

void ElementStore::packData(char* &p) const
{
	for (size_t i = 0; i < data.size(); i++) {
		if (data[i]->names.size()) {
			Esutils::pack(*data[i]->data, p);
		}
	}
}

void ElementStore::unpackData(const char* &p)
{
	for (size_t i = 0; i < data.size(); i++) {
		Esutils::unpack(*data[i]->data, p);
	}
}

ElementStore::~ElementStore()
{
	if (IDs == NULL) { delete IDs; }
	if (nodes == NULL) { delete nodes; }

	if (body == NULL) { delete body; }
	if (material == NULL) { delete material; }
	if (regions == NULL) { delete regions; }
	if (epointers == NULL) { delete epointers; }

	if (neighbors == NULL) { delete neighbors; }
}

void ElementStore::store(const std::string &file)
{
	std::ofstream os(file + std::to_string(environment->MPIrank) + ".txt");

	Store::storedata(os, "IDs", IDs);

	Store::storedata(os, "nodes", nodes);

	Store::storedata(os, "body", body);
	Store::storedata(os, "material", material);
	Store::storedata(os, "epointers", epointers);
}

void ElementStore::permute(const std::vector<eslocal> &permutation, const std::vector<size_t> &distribution)
{
	this->distribution = distribution;

	if (IDs != NULL) { IDs->permute(permutation, distribution); }

	if (nodes != NULL) { nodes->permute(permutation, distribution); }

	if (body != NULL) { body->permute(permutation, distribution); }
	if (material != NULL) { material->permute(permutation, distribution); }
	if (regions != NULL) { regions->permute(permutation, distribution); }

	if (epointers != NULL) {
		size_t threads = environment->OMP_NUM_THREADS;
		if (threads > 1) {
			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				for (size_t i = epointers->datatarray().distribution()[t]; i < epointers->datatarray().distribution()[t + 1]; ++i) {
					epointers->datatarray()[i] = _eclasses[0] + (epointers->datatarray()[i] - _eclasses[t]);
				}
			}

			epointers->permute(permutation, distribution);

			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				for (size_t i = this->distribution[t]; i < this->distribution[t + 1]; ++i) {
					epointers->datatarray()[i] = _eclasses[t] + (epointers->datatarray()[i] - _eclasses[0]);
				}
			}
		} else {
			epointers->permute(permutation, distribution);
		}
	}

	if (neighbors != NULL) { neighbors->permute(permutation, distribution); }
}

ElementData* ElementStore::appendData(int dimension, const std::vector<std::string> &names)
{
	data.push_back(new ElementData(dimension, names));
	if (names.size() < 2) {
		data.back()->data->resize(size);
	} else {
		data.back()->data->resize((names.size() - 1) * size);
	}
	return data.back();
}

ElementData* ElementStore::appendData(int dimension, const std::vector<std::string> &names, std::vector<double> &data)
{
	this->data.push_back(new ElementData(dimension, names, &data));
	return this->data.back();
}

std::vector<eslocal> ElementStore::gatherElementsProcDistribution()
{
	return Store::gatherDistribution(size);
}

std::vector<eslocal> ElementStore::gatherDomainsProcDistribution()
{
	return Store::gatherDistribution(ndomains);
}

std::vector<eslocal> ElementStore::gatherDomainsDistribution()
{
	return Store::gatherDistribution(domainDistribution, firstDomain);
}

std::vector<eslocal> ElementStore::gatherElementsDistribution()
{
	return Store::gatherDistribution(elementsDistribution, IDs->datatarray().front());
}

std::vector<eslocal> ElementStore::gatherClustersDistribution()
{
	return Store::gatherDistribution(nclusters);
}


ElementData::ElementData(int dimension)
: dimension(dimension), _delete(true)
{
	data = new std::vector<double>();
}

ElementData::ElementData(int dimension, const std::vector<std::string> &names, std::vector<double> *data)
: dimension(dimension), names(names), data(data), _delete(false)
{
	if (data == NULL) {
		this->data = new std::vector<double>();
	}
}

ElementData::ElementData(ElementData &&other)
: dimension(std::move(dimension)), names(std::move(other.names)), data(other.data),
  _delete(std::move(other._delete))
{
	other.data = NULL;
	other._delete = false;
}

ElementData::ElementData(const ElementData &other)
: dimension(std::move(dimension)), names(other.names), data(other.data),
  _delete(other._delete)
{
	if (_delete) {
		data = new std::vector<double>(*other.data);
	}
}

ElementData& ElementData::operator=(const ElementData &other)
{
	if (this != &other) {
		dimension = other.dimension;
		names = other.names;
		data = other.data;
		_delete = other._delete;
		if (_delete) {
			data = new std::vector<double>(*other.data);
		}
	}
	return *this;
}

ElementData::~ElementData()
{
	if (_delete && data != NULL) {
		delete data;
	}
}
