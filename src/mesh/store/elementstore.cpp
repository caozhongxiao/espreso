
#include "store.h"
#include "elementstore.h"
#include "statisticsstore.h"

#include "esinfo/mpiinfo.h"
#include "esinfo/envinfo.h"

#include "mesh/mesh.h"
#include "mesh/elements/element.h"

#include "basis/containers/point.h"
#include "basis/containers/serializededata.h"


using namespace espreso;

ElementStore::ElementStore(std::vector<Element*> &eclasses)
: dimension(3),
  size(0),
  distribution({0, 0}),

  IDs(NULL),
  procNodes(NULL),
  domainNodes(NULL),
  centers(NULL),

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
	if (procNodes == NULL || epointers == NULL) {
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
			procNodes->packedSize() +
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
	procNodes->pack(p);
	if (epointers != NULL) {
		std::vector<int> eindices;
		eindices.reserve(epointers->datatarray().size());

		size_t threads = info::env::OMP_NUM_THREADS;
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
	if (procNodes == NULL) {
		procNodes = new serializededata<esint, esint>(tarray<esint>(1, 0), tarray<esint>(1, 0));
	}
	if (epointers == NULL) {
		epointers = new serializededata<esint, Element*>(1, tarray<Element*>(1, 0));
	}

	Esutils::unpack(size, p);
	procNodes->unpack(p);
	if (epointers != NULL) {
		std::vector<int> eindices;
		Esutils::unpack(eindices, p);
		if (epointers != NULL) {
			delete epointers;
		}
		epointers = new serializededata<esint, Element*>(1, tarray<Element*>(1, size));
		for (esint i = 0; i < size; ++i) {
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
			data.push_back(new ElementData(dimension, {}));
		}
		Esutils::unpack(data[i]->names, p);
	}
}

size_t ElementStore::packedDataSize() const
{
	size_t size = 1;
	for (size_t i = 0; i < data.size(); i++) {
		if (data[i]->names.size()) {
			size += Esutils::packedSize(data[i]->data);
		}
	}
	return size;
}

void ElementStore::packData(char* &p) const
{
	for (size_t i = 0; i < data.size(); i++) {
		if (data[i]->names.size()) {
			Esutils::pack(data[i]->data, p);
		}
	}
}

void ElementStore::unpackData(const char* &p)
{
	for (size_t i = 0; i < data.size(); i++) {
		Esutils::unpack(data[i]->data, p);
	}
}

ElementStore::~ElementStore()
{
	if (IDs == NULL) { delete IDs; }
	if (procNodes == NULL) { delete procNodes; }
	if (domainNodes == NULL) { delete domainNodes; }
	if (centers == NULL) { delete centers; }

	if (body == NULL) { delete body; }
	if (material == NULL) { delete material; }
	if (regions == NULL) { delete regions; }
	if (epointers == NULL) { delete epointers; }

	if (neighbors == NULL) { delete neighbors; }
}

void ElementStore::store(const std::string &file)
{
	std::ofstream os(file + std::to_string(info::mpi::rank) + ".txt");

	Store::storedata(os, "IDs", IDs);
	Store::storedata(os, "nodes", procNodes);
	Store::storedata(os, "nodes", domainNodes);
	Store::storedata(os, "centers", centers);

	Store::storedata(os, "body", body);
	Store::storedata(os, "material", material);
	Store::storedata(os, "regions", regions);
	Store::storedata(os, "epointers", epointers);

	Store::storedata(os, "neighbors", neighbors);
}

void ElementStore::permute(const std::vector<esint> &permutation, const std::vector<size_t> &distribution)
{
	this->distribution = distribution;

	if (IDs != NULL) { IDs->permute(permutation, distribution); }
	if (procNodes != NULL) { procNodes->permute(permutation, distribution); }
	if (domainNodes != NULL) { domainNodes->permute(permutation, distribution); }
	if (centers != NULL) { centers->permute(permutation, distribution); }

	if (body != NULL) { body->permute(permutation, distribution); }
	if (material != NULL) { material->permute(permutation, distribution); }
	if (regions != NULL) { regions->permute(permutation, distribution); }

	if (epointers != NULL) {
		size_t threads = info::env::OMP_NUM_THREADS;
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

	// TODO: permute data
}

void ElementStore::reindex(const serializededata<esint, esint> *nIDs)
{
	size_t threads = info::env::OMP_NUM_THREADS;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (auto n = procNodes->begin(t)->begin(); n != procNodes->end(t)->begin(); ++n) {
			*n = std::lower_bound(nIDs->datatarray().begin(), nIDs->datatarray().end(), *n) - nIDs->datatarray().begin();
		}
	}
}

ElementData* ElementStore::appendData(int dimension, const std::vector<std::string> &names)
{
	this->data.push_back(new ElementData(dimension, names));
	data.back()->data.resize(size * dimension);
	return this->data.back();
}

std::vector<esint> ElementStore::gatherElementsProcDistribution()
{
	return Store::gatherDistribution(size);
}

std::vector<esint> ElementStore::gatherDomainsProcDistribution()
{
	return Store::gatherDistribution(ndomains);
}

std::vector<esint> ElementStore::gatherDomainsDistribution()
{
	return Store::gatherDistribution(domainDistribution, firstDomain);
}

std::vector<esint> ElementStore::gatherElementsDistribution()
{
	return Store::gatherDistribution(elementsDistribution, IDs->datatarray().front());
}

std::vector<esint> ElementStore::gatherClustersDistribution()
{
	return Store::gatherDistribution(nclusters);
}

ElementData::ElementData(int dimension, const std::vector<std::string> &names)
: dimension(dimension), names(names)
{

}

void ElementData::statistics(const tarray<esint> &elements, esint totalsize, Statistics *statistics)
{
	for (size_t d = 0; d <= names.size(); d++) {
		(statistics + d)->reset();
	}

	if (names.size() == 1) {
		for (auto e = elements.begin(); e != elements.end(); ++e) {
			statistics->min = std::min(statistics->min, data[*e * dimension]);
			statistics->max = std::max(statistics->max, data[*e * dimension]);
			statistics->avg += data[*e * dimension];
			statistics->norm += data[*e * dimension] * data[*e * dimension];
		}
	} else {
		for (auto e = elements.begin(); e != elements.end(); ++e) {
			double value = 0;
			for (int d = 0; d < dimension; d++) {
				value += data[*e * dimension + d] * data[*e * dimension + d];
				(statistics + d + 1)->min = std::min((statistics + d + 1)->min, data[*e * dimension + d]);
				(statistics + d + 1)->max = std::max((statistics + d + 1)->max, data[*e * dimension + d]);
				(statistics + d + 1)->avg += data[*e * dimension + d];
				(statistics + d + 1)->norm += data[*e * dimension + d] * data[*e * dimension + d];
			}
			value = std::sqrt(value);
			statistics->min = std::min(statistics->min, value);
			statistics->max = std::max(statistics->max, value);
			statistics->avg += value;
			statistics->norm += value * value;
		}
	}

	std::vector<Statistics> global(names.size());
	MPI_Allreduce(statistics, global.data(), sizeof(Statistics) * names.size(), MPI_BYTE, MPITools::operations().mergeStatistics, info::mpi::comm);
	memcpy(statistics, global.data(), sizeof(Statistics) * names.size());

	for (size_t i = 0; i < names.size(); i++) {
		(statistics + i)->avg /= totalsize;
		(statistics + i)->norm = std::sqrt((statistics + i)->norm);
	}
}














