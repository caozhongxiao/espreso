
#ifndef SRC_MESH_STORE_NODESTORE_H_
#define SRC_MESH_STORE_NODESTORE_H_

#include <cstddef>
#include <string>
#include <vector>
#include <iostream>

#include "../intervals/processinterval.h"
#include "../intervals/domaininterval.h"

#include "../../basis/containers/point.h"

namespace espreso {

template <typename TData> class tarray;
template <typename TEBoundaries, typename TEData> class serializededata;
struct Statistics;

struct NodeStore;

struct NodeData {
	int dimension;
	std::vector<std::string> names;

	std::vector<double> data;

	NodeData(int dimension, const std::vector<std::string> &names);

	void statistics(const tarray<eslocal> &nodes, eslocal totalsize, Statistics *statistics) const;
	double norm() const;
};

struct NodeStore {
	friend class Mesh;

	void store(const std::string &file);

	void permute(const std::vector<eslocal> &permutation) { permute(permutation, distribution); }
	void permute(const std::vector<eslocal> &permutation, const std::vector<size_t> &distribution);

	NodeData* appendData(int dimension, const std::vector<std::string> &names = {});

	std::vector<eslocal> gatherNodeDistribution();
	std::vector<eslocal> gatherUniqueNodeDistribution();

	eslocal size;
	eslocal uniqueOffset;
	eslocal uniqueSize;
	eslocal uniqueTotalSize;
	std::vector<size_t> distribution;

	serializededata<eslocal, eslocal>* IDs;
	serializededata<eslocal, eslocal>* elements;

	serializededata<eslocal, Point>* originCoordinates;
	serializededata<eslocal, Point>* coordinates;
	serializededata<eslocal, int>* ranks;

	std::vector<eslocal> externalIntervals;
	std::vector<ProcessInterval> pintervals;
	std::vector<std::vector<DomainInterval> > dintervals;
	serializededata<eslocal, eslocal>* idomains;
	serializededata<eslocal, int>* iranks;

	Point center;
	std::vector<Point> dcenter;
	Point min, max, lmin, lmax;

	std::vector<NodeData*> data;

	size_t packedSize() const;
	void pack(char* &p) const;
	void unpack(const char* &p);

	size_t packedDataSize() const;
	void packData(char* &p) const;
	void unpackData(const char* &p);

	NodeStore();
	~NodeStore();
};

}


#endif /* SRC_MESH_STORE_NODESTORE_H_ */
