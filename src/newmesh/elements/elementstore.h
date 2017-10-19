
#ifndef SRC_NEWMESH_ELEMENTS_ELEMENTSTORE_H_
#define SRC_NEWMESH_ELEMENTS_ELEMENTSTORE_H_

#include <cstddef>
#include <string>
#include <vector>

namespace espreso {

template <typename TEBoundaries, typename TEData> class serializededata;
struct Point;
struct NewElement;

struct ElementStore {

	void store(const std::string &file);

	void sort();
	void permute(const std::vector<eslocal> &permutation, const std::vector<size_t> *distribution = NULL);

	std::vector<esglobal> gatherSizes();

	size_t size;
	std::vector<size_t> distribution;

	serializededata<eslocal, esglobal>* IDs;

	serializededata<eslocal, eslocal>* elems;
	serializededata<eslocal, eslocal>* faces;
	serializededata<eslocal, eslocal>* edges;
	serializededata<eslocal, eslocal>* nodes;

	serializededata<eslocal, Point>* coordinates;
	serializededata<eslocal, int>* body;
	serializededata<eslocal, int>* material;
	serializededata<eslocal, NewElement*>* epointers;
	serializededata<eslocal, eslocal>* domains;
	serializededata<eslocal, int>* ranks;

	serializededata<eslocal, esglobal>* dual;
	serializededata<eslocal, eslocal>* decomposedDual;

	ElementStore();
	~ElementStore();
};

}



#endif /* SRC_NEWMESH_ELEMENTS_ELEMENTSTORE_H_ */
