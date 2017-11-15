
#ifndef SRC_BASIS_CONTAINERS_EDATA_H_
#define SRC_BASIS_CONTAINERS_EDATA_H_

#include <cstddef>

namespace espreso {

template <typename TData>
class edata {

	template <typename TEBoundaries, typename TEData> friend class serializededata;
public:
	size_t size() const          { return this->_end - this->_begin; }

	TData& operator[] (size_t n) { return *(this->_begin + n); }
	TData& at(size_t n)          { return *(this->_begin + n); }
	TData* data()                { return   this->_begin; }
	TData& back()                { return *(this->_end - 1); }
	TData& front()               { return * this->_begin; }
	TData* begin()               { return   this->_begin; }
	TData* end()                 { return   this->_end; }

private:
	edata(TData *begin, TData *end): _begin(begin), _end(end) {}
	TData *_begin;
	TData *_end;
};

}



#endif /* SRC_BASIS_CONTAINERS_EDATA_H_ */
