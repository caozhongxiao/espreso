
#ifndef SRC_MESH_STORE_ELEMENTSTORE_H_
#define SRC_MESH_STORE_ELEMENTSTORE_H_

#include <cstddef>
#include <string>
#include <vector>

#include "mesh/intervals/elementsinterval.h"

namespace espreso {

struct Point;
template <typename TData> class tarray;
template <typename TEBoundaries, typename TEData> class serializededata;
struct Statistics;
struct Element;

struct ElementData {
	int dimension;
	std::vector<std::string> names;

	std::vector<double> data;

	ElementData(int dimension, const std::vector<std::string> &names);

	void statistics(const tarray<esint> &elements, esint totalsize, Statistics *statistics);
};

struct ElementStore {

	void store(const std::string &file);

	void permute(const std::vector<esint> &permutation) { permute(permutation, distribution); }
	void permute(const std::vector<esint> &permutation, const std::vector<size_t> &distribution);

	void reindex(const serializededata<esint, esint> *nIDs);

	ElementData* appendData(int dimension, const std::vector<std::string> &names = {});

	std::vector<esint> gatherElementsDistribution();
	std::vector<esint> gatherElementsProcDistribution();

	std::vector<esint> gatherDomainsDistribution();
	std::vector<esint> gatherDomainsProcDistribution();

	std::vector<esint> gatherClustersDistribution();

	esint dimension;
	esint size;
	std::vector<size_t> distribution;

	serializededata<esint, esint>* IDs;
	serializededata<esint, esint>* procNodes;
	serializededata<esint, esint>* domainNodes;
	serializededata<esint, float>* centers;

	serializededata<esint, int>* body;
	serializededata<esint, int>* material;
	serializededata<esint, esint>* regions;
	serializededata<esint, Element*>* epointers;


	serializededata<esint, esint>* neighbors;

	esint firstDomain;
	esint ndomains;
	std::vector<esint> domainDistribution;
	std::vector<esint> elementsDistribution;
	std::vector<int> clusters;
	esint nclusters;

	int regionMaskSize;
	std::vector<esint> ecounters;
	std::vector<ElementsInterval> eintervals;
	std::vector<esint> eintervalsDistribution;

	std::vector<ElementData*> data;

	size_t packedSize() const;
	void pack(char* &p) const;
	void unpack(const char* &p);

	size_t packedDataSize() const;
	void packData(char* &p) const;
	void unpackData(const char* &p);

	ElementStore(std::vector<Element*> &eclasses);
	~ElementStore();

private:
	std::vector<Element*> &_eclasses;
};

}



#endif /* SRC_MESH_STORE_ELEMENTSTORE_H_ */
