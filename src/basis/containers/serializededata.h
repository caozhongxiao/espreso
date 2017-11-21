
#ifndef SRC_BASIS_CONTAINERS_SERIALIZEDEDATA_H_
#define SRC_BASIS_CONTAINERS_SERIALIZEDEDATA_H_

#include <vector>
#include <ostream>

#include "edata.h"
#include "tarray.h"

namespace espreso {

template <typename TEBoundaries, typename TEData>
class serializededata {

public:
	static void balance(size_t esize, std::vector<std::vector<TEData> > &data, const std::vector<size_t> *distribution = NULL)
	{
		std::vector<size_t> _distribution;
		if (distribution == NULL) {
			size_t size = 0;
			for (size_t t = 0; t < data.size(); t++) {
				size += esize * data[t].size();
			}

			_distribution = tarray<TEBoundaries>::distribute(data.size(), size);
		} else {
			_distribution = *distribution;
		}

		for (size_t t = 0, tt = 0; t < data.size(); tt = ++t) {
			while (++tt && data[t].size() < _distribution[t + 1] - _distribution[t]) {
				size_t diff = _distribution[t + 1] - _distribution[t] - data[t].size();
				if (diff < data[tt].size()) {
					data[t].insert(data[t].end(), data[tt].begin(), data[tt].begin() + diff);
					data[tt].erase(data[tt].begin(), data[tt].begin() + diff);
				} else {
					data[t].insert(data[t].end(), data[tt].begin(), data[tt].end());
					data[tt].clear();

				}
			}
			if (data[t].size() > _distribution[t + 1] - _distribution[t]) {
				size_t diff = data[t].size() - (_distribution[t + 1] - _distribution[t]);
				data[t + 1].insert(data[t + 1].begin(), data[t].end() - diff, data[t].end());
				data[t].erase(data[t].end() - diff, data[t].end());
			}
		}
	}

	static void balance(std::vector<std::vector<TEBoundaries> > &boundaries, std::vector<std::vector<TEData> > &data, const std::vector<size_t> *distribution = NULL)
	{
		size_t size = 0;
		std::vector<size_t> sizes(boundaries.size());
		for (size_t t = 0; t < boundaries.size(); t++) {
			sizes[t] = boundaries[t].size();
			size += boundaries[t].size();
		}
		--sizes[0];

		std::vector<size_t> _distribution;
		if (distribution == NULL) {
			_distribution = tarray<TEBoundaries>::distribute(data.size(), size - 1);
		} else {
			_distribution = *distribution;
		}

		for (size_t t = 0, tt = 0; t < boundaries.size(); tt = ++t) {
			while (++tt && sizes[t] < _distribution[t + 1] - _distribution[t]) {
				size_t diff = _distribution[t + 1] - _distribution[t] - sizes[t];
				if (diff < sizes[tt]) {
					size_t ttt = 0;
					while (boundaries[tt - ++ttt].size() == 0);
					size_t ediff = *(boundaries[tt].begin() + diff - 1) - boundaries[tt - ttt].back();
					data[t].insert(data[t].end(), data[tt].begin(), data[tt].begin() + ediff);
					data[tt].erase(data[tt].begin(), data[tt].begin() + ediff);

					boundaries[t].insert(boundaries[t].end(), boundaries[tt].begin(), boundaries[tt].begin() + diff);
					boundaries[tt].erase(boundaries[tt].begin(), boundaries[tt].begin() + diff);
					sizes[t]  += diff;
					sizes[tt] -= diff;
				} else {
					sizes[t]  += boundaries[tt].size();
					sizes[tt] -= boundaries[tt].size();

					data[t].insert(data[t].end(), data[tt].begin(), data[tt].end());
					data[tt].clear();

					boundaries[t].insert(boundaries[t].end(), boundaries[tt].begin(), boundaries[tt].end());
					boundaries[tt].clear();
				}
			}
			if (sizes[t] > _distribution[t + 1] - _distribution[t]) {
				size_t diff = sizes[t] - (_distribution[t + 1] - _distribution[t]);
				size_t ediff = boundaries[t].back() - *(boundaries[t].end() - diff - 1);
				data[t + 1].insert(data[t + 1].begin(), data[t].end() - ediff, data[t].end());
				data[t].erase(data[t].end() - ediff, data[t].end());

				boundaries[t + 1].insert(boundaries[t + 1].begin(), boundaries[t].end() - diff, boundaries[t].end());
				boundaries[t].erase(boundaries[t].end() - diff, boundaries[t].end());
				sizes[t]     -= diff;
				sizes[t + 1] += diff;
			}
		}
	}

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
		TIterator  operator++(int) {TIterator tmp(*static_cast<TIterator*>(this)); operator++(); return tmp; }
		TIterator  operator--(int) {TIterator tmp(*static_cast<TIterator*>(this)); operator--(); return tmp; }
		template <typename TType> TIterator  operator+ (TType n) { return TIterator(*static_cast<TIterator*>(this)).move( n); }
		template <typename TType> TIterator  operator- (TType n) { return TIterator(*static_cast<TIterator*>(this)).move(-n); }
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
			const_iterator it = _constiterator.front();
			return _edata.size() / it->size();
		}
	}

	void permute(const std::vector<eslocal> &permutation, const std::vector<size_t> &distribution)
	{
		if (_eboundaries.size()) {
			permuteNonUniformData(permutation, distribution);
		} else {
			permuteUniformData(permutation, distribution);
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

	tarray<TEBoundaries>&       boundarytaaray()       { return _eboundaries; }
	const tarray<TEBoundaries>& boundarytaaray() const { return _eboundaries; }
	tarray<TEData>&             datatarray()           { return _edata; }
	const tarray<TEData>&       datatarray()     const { return _edata; }

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

	void permuteUniformData(const std::vector<eslocal> &permutation, const std::vector<size_t> &distribution)
	{
		std::vector<std::vector<TEData> > pdata(threads());

		#pragma omp parallel for
		for (size_t t = 0; t < threads(); t++) {
			for (size_t i = distribution[t]; i < distribution[t + 1]; ++i) {
				pdata[t].push_back(_edata.data()[permutation[i]]);
			}
		}

		_edata = tarray<TEData>(pdata);
		inititerators(_iterator.front()->end() - _iterator.front()->begin());
	}

	void permuteNonUniformData(const std::vector<eslocal> &permutation, const std::vector<size_t> &providedDistribution)
	{
		std::vector<std::vector<TEBoundaries> > pboundaries(threads());
		std::vector<std::vector<TEData> > pdata(threads());
		std::vector<size_t> distribution = providedDistribution;
		if (_eboundaries.distribution().back() > distribution.back()) {
			for (size_t t = 1; t <= threads(); t++) {
				distribution[t] += 1;
			}
		}

		pboundaries.front().push_back(0);
		#pragma omp parallel for
		for (size_t t = 0; t < threads(); t++) {
			for (size_t e = (t == 0 ? 1 : distribution[t]); e < distribution[t + 1]; ++e) {
				pboundaries[t].push_back(_eboundaries.data()[permutation[e - 1] + 1] - _eboundaries.data()[permutation[e - 1]]);
				pdata[t].insert(
						pdata[t].end(),
						_edata.data() + _eboundaries.data()[permutation[e - 1]],
						_edata.data() + _eboundaries.data()[permutation[e - 1] + 1]);
				if (pboundaries[t].size() > 1) {
					pboundaries[t].back() += *(pboundaries[t].end() - 2);
				}
			}
		}

		std::vector<size_t> offsets;
		for (size_t t = 0; t < pboundaries.size(); t++) {
			offsets.push_back(pboundaries[t].size() ? pboundaries[t].back() : 0);
		}

		TEBoundaries sum = 0;
		for (size_t i = 0; i < offsets.size(); i++) {
			TEBoundaries tmp = offsets[i];
			offsets[i] = sum;
			sum += tmp;
		}

		#pragma omp parallel for
		for (size_t t = 0; t < pboundaries.size(); t++) {
			size_t offset = offsets[t];
			for (size_t i = 0; i < pboundaries[t].size(); i++) {
				pboundaries[t][i] += offset;
			}
		}

		_eboundaries = tarray<TEBoundaries>(pboundaries);
		_edata = tarray<TEData>(pdata);
		inititerators();
	}

	tarray<TEBoundaries> _eboundaries;
	tarray<TEData> _edata;

	std::vector<iterator> _iterator;
	std::vector<const_iterator> _constiterator;
};

template <typename TEBoundaries, typename TEData>
struct serializededatainterval {
	serializededatainterval(typename serializededata<TEBoundaries, TEData>::const_iterator begin, typename serializededata<TEBoundaries, TEData>::const_iterator end)
	: _begin(begin), _end(end) {}

	typename serializededata<TEBoundaries, TEData>::const_iterator& begin() { return _begin; }
	typename serializededata<TEBoundaries, TEData>::const_iterator& end() { return _end; }

	typename serializededata<TEBoundaries, TEData>::const_iterator _begin, _end;
};

template <typename TEBoundaries, typename TEData>
std::ostream& operator<< (std::ostream& os, const serializededata<TEBoundaries, TEData> &data)
{
	for(auto e = data.cbegin(); e != data.cend(); ++e) {
		os << *e;
	}
	return os;
}

}


#endif /* SRC_BASIS_CONTAINERS_SERIALIZEDEDATA_H_ */
