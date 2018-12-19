
#ifndef SRC_MESH_STORE_ELEMENTSTORE_H_
#define SRC_MESH_STORE_ELEMENTSTORE_H_

#include <cstddef>
#include <string>
#include <vector>

#include "../intervals/elementsinterval.h"

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

	void statistics(const tarray<eslocal> &elements, eslocal totalsize, Statistics *statistics);
};

struct ElementStore {

	void store(const std::string &file);

	void permute(const std::vector<eslocal> &permutation) { permute(permutation, distribution); }
	void permute(const std::vector<eslocal> &permutation, const std::vector<size_t> &distribution);

	void reindex(const serializededata<eslocal, eslocal> *nIDs);

	ElementData* appendData(int dimension, const std::vector<std::string> &names = {});

	std::vector<eslocal> gatherElementsDistribution();
	std::vector<eslocal> gatherElementsProcDistribution();

	std::vector<eslocal> gatherDomainsDistribution();
	std::vector<eslocal> gatherDomainsProcDistribution();

	std::vector<eslocal> gatherClustersDistribution();

	eslocal dimension;
	eslocal size;
	std::vector<size_t> distribution;

	serializededata<eslocal, eslocal>* IDs;
	serializededata<eslocal, eslocal>* procNodes;
	serializededata<eslocal, eslocal>* domainNodes;
	serializededata<eslocal, double>* centers;

	serializededata<eslocal, int>* body;
	serializededata<eslocal, int>* material;
	serializededata<eslocal, eslocal>* regions;
	serializededata<eslocal, Element*>* epointers;


	serializededata<eslocal, eslocal>* neighbors;

	eslocal firstDomain;
	eslocal ndomains;
	std::vector<eslocal> domainDistribution;
	std::vector<eslocal> elementsDistribution;
	std::vector<int> clusters;
	eslocal nclusters;

	int regionMaskSize;
	std::vector<eslocal> ecounters;
	std::vector<ElementsInterval> eintervals;
	std::vector<eslocal> eintervalsDistribution;

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
