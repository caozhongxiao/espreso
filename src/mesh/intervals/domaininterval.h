
#ifndef SRC_MESH_INTERVALS_DOMAININTERVAL_H_
#define SRC_MESH_INTERVALS_DOMAININTERVAL_H_

namespace espreso {

struct DomainInterval {
	esint begin, end;
	esint pindex, DOFOffset;

	DomainInterval(esint begin, esint end, esint pindex, esint DOFOffset): begin(begin), end(end), pindex(pindex), DOFOffset(DOFOffset) {}
};

}


#endif /* SRC_MESH_INTERVALS_DOMAININTERVAL_H_ */
