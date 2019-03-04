
#ifndef SRC_BASIS_UTILITIES_DEBUGPRINT_H_
#define SRC_BASIS_UTILITIES_DEBUGPRINT_H_

#include "basis/containers/point.h"
#include "basis/containers/serializededata.h"
#include "basis/matrices/denseMatrix.h"
#include "solver/generic/SparseMatrix.h"

#include <ostream>
#include <vector>
#include <iomanip>

namespace espreso {

inline std::ostream& operator<<(std::ostream& os, const Point &p) {
	os << "<" << p.x << " " << p.y << " " << p.z << ">\n";
	return os;
}

inline std::ostream& operator<<(std::ostream& os, const Matrix &m)
{
	for (size_t i = 0; i < m.rows(); i++) {
		for (size_t j = 0; j < m.columns(); j++) {
			os << m(i, j) << " ";
		}
		os << std::endl;
	}
	os << std::endl;
	return os;
}

inline std::ostream& operator<<(std::ostream& os, const SparseMatrix &m)
{
	os << std::setw(6) << m.rows << " " << m.cols << " " << m.nnz << "\n";

	SparseMatrix s = m;
	if (s.CSR_J_col_indices.size()) {
		s.ConvertToCOO(1);
	}
	if (s.dense_values.size()) {
		s.ConvertDenseToCSR(0);
		s.ConvertToCOO(1);
	}

	for (esint i = 0; i < s.nnz; i++) {
		os << std::setw(6) << s.I_row_indices[i] << " ";
		os << std::setw(6) << s.J_col_indices[i] << " ";
		os << std::setw(20) << std::scientific << s.V_values[i] << "\n";
	}
	return os;
}

template<typename T1, typename T2>
std::ostream& operator<<(std::ostream& os, const std::pair<T1, T2> &v)
{
	os << "<" << v.first << ":" << v.second << ">\n";
	return os;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T> &v)
{
	for(size_t i = 0; i < v.size(); ++i) {
		os << std::setw(17) << std::scientific << v[i] << "\n";
	}
	os << "\n";
	return os;
}

template <typename TData>
std::ostream& operator<<(std::ostream& os, edata<TData> &data)
{
	os << "[ ";
	for (auto i = data.begin(); i != data.end(); ++i) {
		os << *i << " ";
	}
	os << "]\n";
	return os;
}

template <typename TEBoundaries, typename TEData>
std::ostream& operator<<(std::ostream& os, const serializededata<TEBoundaries, TEData> &data)
{
	size_t i = 0;
	for(auto e = data.cbegin(); e != data.cend(); ++e, ++i) {
		os << i << ": " << *e;
	}
	return os;
}

}



#endif /* SRC_BASIS_UTILITIES_DEBUGPRINT_H_ */
