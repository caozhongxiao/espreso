
#ifndef SRC_MESH_INTERVALS_ELEMENTSINTERVAL_H_
#define SRC_MESH_INTERVALS_ELEMENTSINTERVAL_H_

namespace espreso {

struct ElementsInterval {
	eslocal begin, end;
	int code;

	ElementsInterval(): begin(0), end(0), code(-1) {}
	ElementsInterval(eslocal begin, eslocal end): begin(begin), end(end), code(0) {}
	ElementsInterval(eslocal begin, eslocal end, int code) : begin(begin), end(end), code(code) {}

	bool operator==(const ElementsInterval &other) const { return begin == other.begin && end == other.end && code == other.code; }
};

}


#endif /* SRC_MESH_INTERVALS_ELEMENTSINTERVAL_H_ */
