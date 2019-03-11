
#ifndef SRC_BASIS_UTILITIES_PRINT_H_
#define SRC_BASIS_UTILITIES_PRINT_H_

#include "basis/containers/point.h"
#include "basis/containers/serializededata.h"
#include "basis/matrices/denseMatrix.h"
#include "solver/generic/SparseMatrix.h"

#include <iostream>
#include <ostream>
#include <vector>

namespace espreso {

inline std::ostream& operator<<(std::ostream& os, const Point &p) {
	os << p.x << " " << p.y << " " << p.z;
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
	os << m.rows << " " << m.cols << " " << m.nnz << "\n";

	os.precision(15);

	SparseMatrix s = m;
	if (s.CSR_J_col_indices.size()) {
		s.ConvertToCOO(1);
	}
	if (s.dense_values.size()) {
		s.ConvertDenseToCSR(0);
		s.ConvertToCOO(1);
	}

	for (esint i = 0; i < s.nnz; i++) {
		os << s.I_row_indices[i] << " ";
		os << s.J_col_indices[i] << " ";
		os << std::scientific << s.V_values[i] << "\n";
	}
	return os;
}

template<typename T1, typename T2>
std::ostream& operator<< (std::ostream& os, const std::pair<T1, T2> &v)
{
	os << "<" << v.first << ":" << v.second << ">";
	return os;
}

template<typename T>
std::ostream& operator<< (std::ostream& os, const std::vector<T> &v)
{
	for(size_t i = 0; i < v.size(); ++i) {
		os << v[i] << " ";
	}
	os << "\n";
	return os;
}

template <typename TData>
std::ostream& operator<< (std::ostream& os, edata<TData> &data)
{
	os << "[ ";
	for (auto i = data.begin(); i != data.end(); ++i) {
		os << *i << " ";
	}
	os << "]";
	return os;
}

template <typename TEBoundaries, typename TEData>
std::ostream& operator<< (std::ostream& os, const serializededata<TEBoundaries, TEData> &data)
{
	for(auto e = data.cbegin(); e != data.cend(); ++e) {
		os << *e;
	}
	return os;
}

}

#endif /* SRC_BASIS_UTILITIES_PRINT_H_ */
