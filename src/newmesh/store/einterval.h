
#ifndef SRC_NEWMESH_STORE_EINTERVAL_H_
#define SRC_NEWMESH_STORE_EINTERVAL_H_

namespace espreso {

struct EInterval {
	eslocal begin, end;
	std::vector<int> neighbors;

	EInterval(eslocal begin, eslocal end, std::vector<int> &&neighbors)
	: begin(begin), end(end), neighbors(std::move(neighbors)) {};
};

}



#endif /* SRC_NEWMESH_STORE_EINTERVAL_H_ */
