
#ifndef SRC_NEWMESH_ELEMENTS_ELEMENT_H_
#define SRC_NEWMESH_ELEMENTS_ELEMENT_H_

#include <cstddef>

namespace espreso {

template <typename TEData> class edata;
class ElementStore;

class NewElement {

public:
	NewElement(size_t index, ElementStore *store): _index(index), _store(store) {}

	edata<eslocal> nodes();
	edata<const eslocal> nodes() const;

protected:
	size_t _index;
	ElementStore *_store;
};

}


#endif /* SRC_NEWMESH_ELEMENTS_ELEMENT_H_ */
