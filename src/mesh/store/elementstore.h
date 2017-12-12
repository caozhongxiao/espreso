
#ifndef SRC_MESH_STORE_ELEMENTSTORE_H_
#define SRC_MESH_STORE_ELEMENTSTORE_H_

#include <cstddef>
#include <string>
#include <vector>

#include "../intervals/elementsinterval.h"

namespace espreso {

template <typename TEBoundaries, typename TEData> class serializededata;
struct Element;

struct ElementData {
	friend class Mesh;

	std::vector<std::string> names;

	std::vector<double> *data;

	size_t packedSize() const;
	void pack(char* &p) const;
	void unpack(const char* &p);

	ElementData();
	ElementData(const std::vector<std::string> &names, std::vector<double> *data = NULL);
	ElementData(std::vector<double> *data);

	ElementData(ElementData &&other);
	ElementData(const ElementData &other);
	ElementData& operator=(const ElementData&other);
	~ElementData();

protected:
	bool _delete;
};

struct ElementStore {

	void store(const std::string &file);

	void permute(const std::vector<eslocal> &permutation) { permute(permutation, distribution); }
	void permute(const std::vector<eslocal> &permutation, const std::vector<size_t> &distribution);

	ElementData* appendData(const std::vector<std::string> &names);
	ElementData* appendData(const std::vector<std::string> &names, std::vector<double> &data);

	std::vector<eslocal> gatherElementsDistribution();
	std::vector<eslocal> gatherElementsProcDistribution();

	std::vector<eslocal> gatherDomainsDistribution();
	std::vector<eslocal> gatherDomainsProcDistribution();

	std::vector<eslocal> gatherClustersDistribution();

	eslocal size;
	std::vector<size_t> distribution;

	serializededata<eslocal, eslocal>* IDs;
	serializededata<eslocal, eslocal>* nodes;

	serializededata<eslocal, int>* body;
	serializededata<eslocal, int>* material;
	serializededata<eslocal, int>* regions;
	serializededata<eslocal, Element*>* epointers;

	serializededata<eslocal, eslocal>* dual;
	serializededata<eslocal, eslocal>* decomposedDual;

	eslocal firstDomain;
	eslocal ndomains;
	std::vector<eslocal> domainDistribution;
	std::vector<eslocal> elementsDistribution;
	std::vector<int> clusters;
	eslocal nclusters;

	int regionMaskSize;
	std::vector<eslocal> ecounters;
	std::vector<ElementsInterval> eintervals;

	std::vector<ElementData*> data;

	size_t packedSize() const;
	void pack(char* &p) const;
	void unpack(const char* &p);

	ElementStore(std::vector<Element*> &eclasses);
	~ElementStore();

private:
	std::vector<Element*> &_eclasses;
};

}



#endif /* SRC_MESH_STORE_ELEMENTSTORE_H_ */
