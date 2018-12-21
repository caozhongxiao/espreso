
#ifndef SRC_MESH_INTERVALS_PROCESSINTERVAL_H_
#define SRC_MESH_INTERVALS_PROCESSINTERVAL_H_

namespace espreso {

struct ProcessInterval {
	esint begin, end;
	esint sourceProcess, globalOffset;

	ProcessInterval(): begin(0), end(0), sourceProcess(0), globalOffset(0) {}
	ProcessInterval(esint begin, esint end): begin(begin), end(end), sourceProcess(0), globalOffset(0) {}
	ProcessInterval(esint begin, esint end, esint sourceProcess, esint globalOffset)
	: begin(begin), end(end), sourceProcess(sourceProcess), globalOffset(globalOffset) {}
};

}


#endif /* SRC_MESH_INTERVALS_PROCESSINTERVAL_H_ */
