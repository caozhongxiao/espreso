
#include "domainscomposer.h"

#include "../../../instance.h"

#include "../../../../solver/generic/SparseMatrix.h"

using namespace espreso;

void DomainsComposer::clearMatrices(Matrices matrices, eslocal domain)
{
	if (matrices & Matrices::K) {
		std::fill(_instance.K[domain].CSR_V_values.begin(), _instance.K[domain].CSR_V_values.end(), 0);
	}
	if (matrices & Matrices::M) {
		std::fill(_instance.M[domain].CSR_V_values.begin(), _instance.M[domain].CSR_V_values.end(), 0);
	}
	if (matrices & Matrices::f) {
		std::fill(_instance.f[domain].begin(), _instance.f[domain].end(), 0);
	}
	if (matrices & Matrices::R) {
		std::fill(_instance.R[domain].begin(), _instance.R[domain].end(), 0);
	}
}


