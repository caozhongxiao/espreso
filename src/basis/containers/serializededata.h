
#ifndef SRC_BASIS_CONTAINERS_SERIALIZEDEDATA_H_
#define SRC_BASIS_CONTAINERS_SERIALIZEDEDATA_H_

#include <vector>

#include "edata.h"
#include "tarray.h"

namespace espreso {

template <typename TEBoundaries, typename TEData>
class serializededata {

private:
	template<class TIterator, typename TIteratorEData>
	class iterator_base {

	public:
		bool operator< (const TIterator &other) const { return _edata._begin <  other._edata._begin; }
		bool operator> (const TIterator &other) const { return _edata._begin >  other._edata._begin; }
		bool operator<=(const TIterator &other) const { return _edata._begin <= other._edata._begin; }
		bool operator>=(const TIterator &other) const { return _edata._begin >= other._edata._begin; }
		bool operator==(const TIterator &other) const { return _edata._begin == other._edata._begin; }
		bool operator!=(const TIterator &other) const { return _edata._begin != other._edata._begin; }

		TIterator& operator++() { return move( 1); }
		TIterator& operator--() { return move(-1); }
		TIterator  operator++(int) {TIterator tmp(static_cast<TIterator*>(this)); operator++(); return tmp; }
		TIterator  operator--(int) {TIterator tmp(static_cast<TIterator*>(this)); operator--(); return tmp; }
		template <typename TType> TIterator& operator+ (TType n) { return move( n); }
		template <typename TType> TIterator& operator- (TType n) { return move(-n); }
		template <typename TType> TIterator& operator+=(TType n) { return move( n); }
		template <typename TType> TIterator& operator-=(TType n) { return move(-n); }

	protected:
		iterator_base(const TEBoundaries *begin, const TEBoundaries *element, const TEBoundaries *end, TIteratorEData *edata)
		: _element(begin), _end(end), _edata(edata, edata)
		{
			move(element - begin);
		}

		iterator_base(size_t edatasize, TIteratorEData *edata)
		: _element(NULL), _end(NULL), _edata(edata, edata)
		{
			_edata._end += edatasize;
		}

		template <typename TType>
		TIterator& move(TType n)
		{
			if (_element == NULL) {
				size_t size = _edata._end - _edata._begin;
				_edata._begin += n * size;
				_edata._end += n * size;
			} else {
				_edata._begin += *(_element + n) - *(_element);
				_element += n;
				if (_element != _end) {
					_edata._end = _edata._begin + *(_element + 1) - *_element;
				} else {
					_edata._end = _edata._begin;
				}
			}
			return static_cast<TIterator&>(*this);
		}

		const TEBoundaries* _element;
		const TEBoundaries* _end;
		edata<TIteratorEData> _edata;
	};

public:

	class iterator: public iterator_base<iterator, TEData>
	{
		friend class serializededata<TEBoundaries, TEData>;
	public:
		edata<TEData>& operator*()  { return  this->_edata; }
		edata<TEData>* operator->() { return &this->_edata; }

	private:
		iterator(TEBoundaries *begin, TEBoundaries *element, TEBoundaries *end, TEData *edata)
		: iterator_base<iterator, TEData>(begin, element, end, edata) { }
		iterator(size_t edatasize, TEData *edata)
		: iterator_base<iterator, TEData>(edatasize, edata) { }
	};

	class const_iterator: public iterator_base<const_iterator, const TEData>
	{
		friend class serializededata<TEBoundaries, TEData>;
	public:
		edata<const TEData>& operator*()  { return  this->_edata; }
		edata<const TEData>* operator->() { return &this->_edata; }

	private:
		const_iterator(const TEBoundaries *begin, const TEBoundaries *element, const TEBoundaries *end, const TEData *edata)
		: iterator_base<const_iterator, const TEData>(begin, element, end, edata) { }
		const_iterator(size_t edatasize, TEData *edata)
		: iterator_base<const_iterator, const TEData>(edatasize, edata) { }
	};

	// data are uniform
	serializededata(size_t edatasize, tarray<TEData> &&edata)
	: _eboundaries(0, 0), _edata(std::move(edata)) { inititerators(edatasize); }

	// data are non-uniform
	serializededata(tarray<TEBoundaries> &&eboundaries, tarray<TEData> &&edata)
	: _eboundaries(std::move(eboundaries)), _edata(std::move(edata)) { inititerators(); }

	serializededata(const serializededata<TEBoundaries, TEData> &other)
	: _eboundaries(other._eboundaries), _edata(other._edata) { inititerators(); }
	serializededata(serializededata<TEBoundaries, TEData> &&other)
	: _eboundaries(std::move(other._eboundaries)), _edata(std::move(other._edata)) { inititerators(); }

	serializededata<TEBoundaries, TEData>& operator=(const serializededata<TEBoundaries, TEData> &other)
	{
		if (this != &other) {
			_eboundaries = other._eboundaries;
			_edata = other._edata;
			inititerators();
		}
		return *this;
	}
	serializededata<TEBoundaries, TEData>& operator=(serializededata<TEBoundaries, TEData> &&other)
	{
		if (this != &other) {
			_eboundaries = std::move(other._eboundaries);
			_edata = std::move(other._edata);
			inititerators();
		}
		return *this;
	}

	size_t threads() const { return _edata.threads(); }
	size_t structures() const
	{
		if (_eboundaries.size()) {
			return _eboundaries.size() - 1;
		} else {
			return _edata.size() / _iterator.front()->size();
		}
	}

	iterator       begin()        { return _iterator.front(); }
	const_iterator begin()  const { return _constiterator.front(); }
	const_iterator cbegin() const { return _constiterator.front(); }
	iterator       end()          { return _iterator.back(); }
	const_iterator end()    const { return _constiterator.back(); }
	const_iterator cend()   const { return _constiterator.back(); }

	iterator       begin (size_t thread)       { return _iterator[thread]; }
	const_iterator begin (size_t thread) const { return _constiterator[thread]; }
	const_iterator cbegin(size_t thread) const { return _constiterator[thread]; }
	iterator       end   (size_t thread)       { return _iterator[thread + 1]; }
	const_iterator end   (size_t thread) const { return _constiterator[thread + 1]; }
	const_iterator cend  (size_t thread) const { return _constiterator[thread + 1]; }

	tarray<TEData>&       data()       { return _edata; }
	const tarray<TEData>& data() const { return _edata; }

private:
	void inititerators()
	{
		_iterator = std::vector<iterator>(threads() + 1, iterator(_eboundaries.begin(), _eboundaries.begin(), _eboundaries.end() - 1, _edata.begin()));
		_constiterator = std::vector<const_iterator>(threads() + 1, const_iterator(_eboundaries.begin(), _eboundaries.begin(), _eboundaries.end() - 1, _edata.begin()));

		for (size_t t = 1; t <= threads(); t++) {
			_iterator[t] += _eboundaries.distribution()[t] - 1;
			_constiterator[t] += _eboundaries.distribution()[t] - 1;
		}
	}

	void inititerators(size_t edatasize)
	{
		_iterator = std::vector<iterator>(threads() + 1, iterator(edatasize, _edata.begin()));
		_constiterator = std::vector<const_iterator>(threads() + 1, const_iterator(edatasize, _edata.begin()));

		for (size_t t = 1; t <= threads(); t++) {
			_iterator[t] += _edata.distribution()[t] / edatasize;
			_constiterator[t] += _edata.distribution()[t] / edatasize;
		}
	}

	tarray<TEBoundaries> _eboundaries;
	tarray<TEData> _edata;

	std::vector<iterator> _iterator;
	std::vector<const_iterator> _constiterator;
};

}


#endif /* SRC_BASIS_CONTAINERS_SERIALIZEDEDATA_H_ */
