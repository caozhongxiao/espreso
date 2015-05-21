#include "esmesh.h"

void test_matrices();
void test_BEM();

int main(int argc, char** argv)
{
	test_matrices();
	return 0;
}

void test_BEM()
{
	int partsCount = 4;
	int fixPointsCount = 4;

	Mesh mesh("matrices/TET/5/elem", "matrices/TET/5/coord", partsCount, fixPointsCount);

	int dimension = mesh.getPartNodesCount(0) * Point::size();

	SparseCSRMatrix K(dimension, dimension);
	SparseCSRMatrix M(dimension, dimension);
	std::vector<double> f(dimension);

	mesh.elasticity(K, M, f, 0);

	mesh.saveVTK("mesh.vtk");

	Boundaries b(mesh);

	BoundaryMesh bMesh;
	mesh.getBoundary(bMesh);

	std::vector<DenseMatrix> K_mat;

	K_mat.reserve(partsCount);
	for (int d = 0; d < partsCount; d++) {
		K_mat.push_back( DenseMatrix (0, 0) );
	}

	for (int d = 0; d < partsCount; d++) {

		bMesh.elasticity(K_mat[d], d);

		std::cout << d << " " << std::endl;
	}

	//bem.saveVTK("bem.vtk");

	//mesh.saveVTK();
}

void fill_matrix(Matrix *m)
{
	m->resize(5, 4);
	m->operator ()(0, 0) = 1;
	m->operator ()(1, 1) = 2;
	m->operator ()(2, 2) = 3;
	m->operator ()(3, 3) = 4;
	m->operator ()(0, 3) = 3.3;
	m->set(4, 3, 5.5);
}

void test_matrices()
{
	std::vector<Matrix*> matrices;

	DenseMatrix d;
	SparseDOKMatrix dok;
	SparseVVPMatrix vvp;

	matrices.push_back(&d);
	fill_matrix(matrices.back());

	matrices.push_back(&dok);
	fill_matrix(matrices.back());

	matrices.push_back(&vvp);
	fill_matrix(matrices.back());
	matrices.back()->set(1, 1, 3);
	matrices.back()->set(1, 1, -3);

	// IJV matrix
	matrices.push_back(new SparseIJVMatrix(d));
	matrices.push_back(new SparseIJVMatrix());
	*dynamic_cast<SparseIJVMatrix*>(matrices.back()) = d;
	matrices.push_back(new SparseIJVMatrix(dok));
	matrices.push_back(new SparseIJVMatrix());
	*dynamic_cast<SparseIJVMatrix*>(matrices.back()) = dok;
	matrices.push_back(new SparseIJVMatrix(vvp));
	matrices.push_back(new SparseIJVMatrix());
	*dynamic_cast<SparseIJVMatrix*>(matrices.back()) = vvp;

	// CSR matrix
	matrices.push_back(new SparseCSRMatrix(d));
	matrices.push_back(new SparseCSRMatrix());
	*dynamic_cast<SparseCSRMatrix*>(matrices.back()) = d;
	matrices.push_back(new SparseCSRMatrix(dok));
	matrices.push_back(new SparseCSRMatrix());
	*dynamic_cast<SparseCSRMatrix*>(matrices.back()) = dok;

	SparseCSRMatrix ccc(d);
	SparseIJVMatrix iii(d);

	matrices.push_back(new DenseMatrix(ccc));
	matrices.push_back(new DenseMatrix(iii));

	vvp.shrink();
	for (size_t i = 1; i < matrices.size(); i++) {
		for (size_t r = 0; r < matrices[i]->rows(); r++) {
			for (size_t c = 0; c < matrices[i]->columns(); c++) {
				const Matrix *m = matrices[i];
				if (matrices[0]->get(r, c) != m->operator ()(r, c)) {
					std::cerr << *matrices[0];
					std::cerr << *matrices[i];
					std::cerr << "Matrix: " << i << ", ";
					std::cerr << "row: " << r << ", column: " << c << " -> ";
					std::cerr << matrices[0]->get(r, c) << " != ";
					std::cerr << m->operator ()(r, c) << "\n";
					return;
				}
			}
		}
	}

	for (size_t i = 3; i < matrices.size(); i++) {
		delete matrices[i];
	}

	SparseDOKMatrix dokA(4, 6);
	SparseDOKMatrix dokB(6, 5);
	SparseDOKMatrix dokResult(4, 5);
	int result[] = { 1, 4, 10, 20, 35 };

	for (int i = 0; i < 4; i++) {
		for (int j = i; j < 5; j++) {
			dokResult(i, j) = result[j - i];
		}
	}

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 6 - i; j++) {
			dokA(i, j + i) = j + 1;
		}
	}

	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 5 - i; j++) {
			dokB(i, j + i) = j + 1;
		}
	}

	SparseCSRMatrix A(dokA);
	SparseCSRMatrix B(dokB);
	SparseCSRMatrix C;

	C.multiply(A, B);
	DenseMatrix denseC(C);

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 5; j++) {
			if (denseC(i, j) != dokResult(i, j)) {
				std::cerr << "CSR A * CSR B is incorrect\n";
				exit(EXIT_FAILURE);
			}
		}
	}

	dokA.transpose();
	A = dokA;
	SparseCSRMatrix D;
	D.multiply(A, B, true);
	DenseMatrix denseD(D);

	std::cout << dokA;
	std::cout << A;

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 5; j++) {
			if (denseD(i, j) != dokResult(i, j)) {
				std::cerr << "trans CSR A * CSR B is incorrect\n";
				exit(EXIT_FAILURE);
			}
		}
	}

}



