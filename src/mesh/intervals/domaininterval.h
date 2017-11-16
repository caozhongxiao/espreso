
#ifndef SRC_MESH_INTERVALS_DOMAININTERVAL_H_
#define SRC_MESH_INTERVALS_DOMAININTERVAL_H_

namespace espreso {

struct DomainInterval {
	eslocal begin, end;
	eslocal pindex;

	DomainInterval(eslocal begin, eslocal end, eslocal pindex): begin(begin), end(end), pindex(pindex) {}
};

}


#endif /* SRC_MESH_INTERVALS_DOMAININTERVAL_H_ */
