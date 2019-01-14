
#ifndef SRC_MESH_STORE_NODESTORE_H_
#define SRC_MESH_STORE_NODESTORE_H_

#include <cstddef>
#include <string>
#include <vector>
#include <iostream>

#include "mesh/intervals/processinterval.h"
#include "mesh/intervals/domaininterval.h"

#include "basis/containers/point.h"

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

	void statistics(const tarray<esint> &nodes, esint totalsize, Statistics *statistics) const;
	double norm() const;
};

struct NodeStore {
	friend class Mesh;

	void store(const std::string &file);

	void permute(const std::vector<esint> &permutation) { permute(permutation, distribution); }
	void permute(const std::vector<esint> &permutation, const std::vector<size_t> &distribution);

	NodeData* appendData(int dimension, const std::vector<std::string> &names = {});

	std::vector<esint> gatherNodeDistribution();
	std::vector<esint> gatherUniqueNodeDistribution();

	esint size;
	esint uniqueOffset;
	esint uniqueSize;
	esint uniqueTotalSize;
	std::vector<size_t> distribution;

	serializededata<esint, esint>* IDs;
	serializededata<esint, esint>* elements;

	serializededata<esint, Point>* originCoordinates;
	serializededata<esint, Point>* coordinates;
	serializededata<esint, int>* ranks;

	std::vector<esint> externalIntervals;
	std::vector<ProcessInterval> pintervals;
	std::vector<std::vector<DomainInterval> > dintervals;
	serializededata<esint, esint>* idomains;
	serializededata<esint, int>* iranks;

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
