
#ifndef SRC_MESH_INTERVALS_PROCESSINTERVAL_H_
#define SRC_MESH_INTERVALS_PROCESSINTERVAL_H_

namespace espreso {

struct ProcessInterval {
	eslocal begin, end;

	ProcessInterval(): begin(0), end(0) {}
	ProcessInterval(eslocal begin, eslocal end): begin(begin), end(end) {}
};

}


#endif /* SRC_MESH_INTERVALS_PROCESSINTERVAL_H_ */
