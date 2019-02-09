
#ifndef SRC_MESH_INTERVALS_DOMAININTERVAL_H_
#define SRC_MESH_INTERVALS_DOMAININTERVAL_H_

namespace espreso {

struct DomainInterval {
	esint begin, end;
	esint pindex, DOFOffset;

	DomainInterval(): begin(0), end(0), pindex(0), DOFOffset() {} // allow unpack
	DomainInterval(esint begin, esint end, esint pindex, esint DOFOffset): begin(begin), end(end), pindex(pindex), DOFOffset(DOFOffset) {}
};

}


#endif /* SRC_MESH_INTERVALS_DOMAININTERVAL_H_ */
