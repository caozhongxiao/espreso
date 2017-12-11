
#ifndef SRC_MESH_STORE_NODESTORE_H_
#define SRC_MESH_STORE_NODESTORE_H_

#include <cstddef>
#include <string>
#include <vector>
#include <iostream>

#include "../intervals/processinterval.h"
#include "../intervals/domaininterval.h"
#include "../intervals/gluinginterval.h"

#include "../../basis/containers/point.h"

namespace espreso {

template <typename TEBoundaries, typename TEData> class serializededata;

struct NodeStore;

struct NodeData {
	friend class Mesh;

	std::vector<std::string> names;

	std::vector<std::vector<double> > *decomposedData;

	std::vector<double> gatheredData;

	size_t packedSize() const;
	void pack(char* &p) const;
	void unpack(const char* &p);

	NodeData();
	NodeData(const std::vector<std::string> &names, std::vector<std::vector<double> > *data = NULL);
	NodeData(std::vector<std::vector<double> > *data);

	NodeData(NodeData &&other);
	NodeData(const NodeData &other);
	NodeData& operator=(const NodeData&other);
	~NodeData();

protected:
	std::vector<std::vector<double> > sBuffer, rBuffer;

	bool _delete;
};

struct TNeighborOffset {
	eslocal process, offset;
};

inline std::ostream& operator<<(std::ostream& os, const TNeighborOffset &neighborOffset)
{
	os << neighborOffset.process << " " << neighborOffset.offset;
	return os;
}

struct NodeStore {
	friend class Mesh;

	void store(const std::string &file);

	void permute(const std::vector<eslocal> &permutation) { permute(permutation, distribution); }
	void permute(const std::vector<eslocal> &permutation, const std::vector<size_t> &distribution);

	NodeData* appendData(const std::vector<std::string> &names);
	NodeData* appendData(const std::vector<std::string> &names, std::vector<std::vector<double> > *data);

	std::vector<eslocal> gatherNodeDistribution();
	std::vector<eslocal> gatherUniqueNodeDistribution();

	eslocal size;
	eslocal uniqueOffset;
	eslocal uniqueSize;
	eslocal uniqueTotalSize;
	std::vector<size_t> distribution;

	serializededata<eslocal, eslocal>* IDs;
	serializededata<eslocal, eslocal>* elements;

	serializededata<eslocal, Point>* coordinates;
	serializededata<eslocal, int>* ranks;

	std::vector<eslocal> externalIntervals;
	std::vector<ProcessInterval> pintervals;
	std::vector<std::vector<DomainInterval> > dintervals;
	std::vector<std::vector<GluingInterval> > gintervals;
	serializededata<eslocal, eslocal>* idomains;
	serializededata<eslocal, TNeighborOffset>* ineighborOffsets;
	serializededata<eslocal, int>* iranks;

	Point center;
	std::vector<Point> dcenter;

	std::vector<NodeData*> data;
	std::vector<eslocal> soffsets, scouters;

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
