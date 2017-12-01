
#ifndef SRC_MESH_INTERVALS_REGIONINTERVAL_H_
#define SRC_MESH_INTERVALS_REGIONINTERVAL_H_

namespace espreso {

struct RegionInterval {
	eslocal begin, end;

	RegionInterval(): begin(0), end(0) {}
	RegionInterval(eslocal begin, eslocal end): begin(begin), end(end) {}
};

}


#endif /* SRC_MESH_INTERVALS_REGIONINTERVAL_H_ */
