
#ifndef SRC_MESH_INTERVALS_PROCESSINTERVAL_H_
#define SRC_MESH_INTERVALS_PROCESSINTERVAL_H_

namespace espreso {

struct ProcessInterval {
	eslocal begin, end;
	eslocal sourceProcess, globalOffset;

	ProcessInterval(): begin(0), end(0), sourceProcess(0), globalOffset(0) {}
	ProcessInterval(eslocal begin, eslocal end): begin(begin), end(end), sourceProcess(0), globalOffset(0) {}
	ProcessInterval(eslocal begin, eslocal end, eslocal sourceProcess, eslocal globalOffset)
	: begin(begin), end(end), sourceProcess(sourceProcess), globalOffset(globalOffset) {}
};

}


#endif /* SRC_MESH_INTERVALS_PROCESSINTERVAL_H_ */
