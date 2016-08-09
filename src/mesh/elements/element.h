#ifndef ELEMENT_H_
#define ELEMENT_H_

#include <vector>
#include <cstring>

#include "esbasis.h"
#include "../settings/setting.h"

namespace espreso {

class Mesh;

class Element
{
	friend class Mesh;

public:
	enum Params {
		MATERIAL,
		CONSTANT,
		COORDINATES,
		BODY,
		PARAMS_SIZE
	};

	inline static bool match(const eslocal *indices, eslocal x, eslocal y)
	{
		return indices[x] == indices[y];
	}

	friend std::ofstream& operator<<(std::ofstream& os, const Element &e);
	friend std::ostream& operator<<(std::ostream& os, const Element &e);

	bool operator<(const Element& other) const
	{
		if (coarseNodes() != other.coarseNodes()) {
			return coarseNodes() < other.coarseNodes();
		}

		std::vector<eslocal> e1(indices(), indices() + coarseNodes());
		std::vector<eslocal> e2(other.indices(), other.indices() + other.coarseNodes());
		std::sort(e1.begin(), e1.end());
		std::sort(e2.begin(), e2.end());
		return e1 < e2;
	}

	bool operator==(const Element& other) const
	{
		if (coarseNodes() != other.coarseNodes()) {
			return false;
		}
		return std::is_permutation(indices(), indices() + coarseNodes(), other.indices());
	}

	virtual ~Element() {};

	virtual const std::vector<DenseMatrix>& dN() const = 0;
	virtual const std::vector<DenseMatrix>& N() const = 0;
	virtual const std::vector<double>& weighFactor() const = 0;

	virtual const std::vector<Property>& elementDOFs() const = 0;
	virtual const std::vector<Property>& faceDOFs() const = 0;
	virtual const std::vector<Property>& edgeDOFs() const = 0;
	virtual const std::vector<Property>& pointDOFs() const = 0;
	virtual const std::vector<Property>& midPointDOFs() const = 0;

	virtual eslocal nCommon() const = 0;
	virtual eslocal vtkCode() const = 0;

	virtual size_t faces() const = 0;
	virtual size_t edges() const = 0;
	virtual size_t nodes() const = 0;
	virtual size_t coarseNodes() const = 0;
	virtual size_t gaussePoints() const = 0;

	virtual Element* face(size_t index) const = 0;
	virtual Element* edge(size_t index) const = 0;
	eslocal& node(size_t index) { return indices()[index]; }
	const eslocal& node(size_t index) const { return indices()[index]; }

	virtual eslocal param(Params param) const =0;
	virtual void param(Params param, eslocal value) =0;

	Settings& settings() { return _settings; }
	const Settings& settings() const { return _settings; }

	std::vector<Evaluator*>& settings(Property property) { return _settings[property]; }
	const std::vector<Evaluator*>& settings(Property property) const { return _settings[property]; }

	std::vector<Element*>& elements() { return _elements; }
	const std::vector<Element*>& elements() const { return _elements; }

	std::vector<eslocal>& domains() { return _domains; }
	const std::vector<eslocal>& domains() const { return _domains; }

	std::vector<eslocal>& clusters() { return _clusters; }
	const std::vector<eslocal>& clusters() const { return _clusters; }

	std::vector<eslocal>& DOFsIndices() { return _DOFsIndices; }
	const std::vector<eslocal>& DOFsIndices() const { return _DOFsIndices; }

	std::vector<eslocal>& neighbourDOFsCounter() { return _neighbourDOFsCounter; }
	const std::vector<eslocal>& neighbourDOFsCounter() const { return _neighbourDOFsCounter; }

	eslocal DOFIndex(eslocal domain, size_t DOFIndex)
	{
		auto it = std::lower_bound(_domains.begin(), _domains.end(), domain);
		auto DOFs = _DOFsIndices.size() / _domains.size();
		return _DOFsIndices[DOFs * (it - _domains.begin()) + DOFIndex];
	}

	size_t numberOfDomainsWithDOF(size_t index)
	{
		size_t n = 0;
		for (size_t d = 0; d < _domains.size(); d++) {
			if (_DOFsIndices[d * _DOFsIndices.size() / _domains.size() + index] != -1) {
				n++;
			}
		}
		return n;
	}

protected:
	virtual eslocal* indices() = 0;
	virtual const eslocal* indices() const = 0;
	virtual std::vector<eslocal> getNeighbours(size_t nodeIndex) const = 0;

	virtual void face(size_t index, Element* face) = 0;
	virtual void edge(size_t index, Element* edge) = 0;

	virtual void fillFaces() = 0;
	virtual void fillEdges() = 0;

	Settings _settings;
	std::vector<Element*> _elements;
	std::vector<eslocal> _domains;
	std::vector<eslocal> _clusters;
	std::vector<eslocal> _DOFsIndices;
	std::vector<eslocal> _neighbourDOFsCounter;
};

inline std::ofstream& espreso::operator<<(std::ofstream& os, const Element &e)
{
	eslocal value = e.vtkCode();
	os.write(reinterpret_cast<const char *>(&value), sizeof(eslocal));
	os.write(reinterpret_cast<const char *>(e.indices()), sizeof(eslocal) * e.nodes());
	return os;
}

inline std::ostream& operator<<(std::ostream& os, const Element &e)
{
	for (size_t i = 0; i < e.nodes(); i++) {
		os << e.node(i) << " ";
	}
	return os;
}

}


#endif /* ELEMENT_H_ */
