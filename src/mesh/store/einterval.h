
#ifndef SRC_MESH_STORE_EINTERVAL_H_
#define SRC_MESH_STORE_EINTERVAL_H_

namespace espreso {

struct EInterval {
	eslocal begin, end;
	eslocal firstDomain, globalDomainOffset, localDomainOffset, ndomains;
	eslocal globalOffset, clusterOffset, domainOffset, LMOffset;
	std::vector<int> neighbors;

	EInterval(eslocal begin, eslocal end, std::vector<int> &&neighbors)
	: begin(begin), end(end),
	  firstDomain(-1), globalDomainOffset(-1), localDomainOffset(-1), ndomains(-1),
	  globalOffset(-1), clusterOffset(-1), domainOffset(-1), LMOffset(-1),
	  neighbors(std::move(neighbors)) {};
};

}



#endif /* SRC_MESH_STORE_EINTERVAL_H_ */
