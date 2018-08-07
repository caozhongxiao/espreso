
#ifndef SRC_MESH_INTERVALS_DOMAININTERVAL_H_
#define SRC_MESH_INTERVALS_DOMAININTERVAL_H_

namespace espreso {

struct DomainInterval {
	eslocal begin, end;
	eslocal pindex, DOFOffset;

	DomainInterval(eslocal begin, eslocal end, eslocal pindex, eslocal DOFOffset): begin(begin), end(end), pindex(pindex), DOFOffset(DOFOffset) {}
};

}


#endif /* SRC_MESH_INTERVALS_DOMAININTERVAL_H_ */
