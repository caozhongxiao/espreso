
#ifndef SRC_NEWMESH_ELEMENTS_ELEMENTSTORE_H_
#define SRC_NEWMESH_ELEMENTS_ELEMENTSTORE_H_

namespace espreso {

template <typename TEBoundaries, typename TEData> class serializededata;

struct ElementStore {

	serializededata<eslocal, eslocal>* indices;
};

}



#endif /* SRC_NEWMESH_ELEMENTS_ELEMENTSTORE_H_ */
