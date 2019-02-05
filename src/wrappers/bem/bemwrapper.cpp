
#include "bemwrapper.h"

#include <cstddef>
#include <cstdio>

#ifdef HAVE_BEM
#include "heatdtn.h"
#endif

namespace espreso {

#ifdef HAVE_BEM
struct BEMData
{
	bem4i::bem4idata<esint, double> * data;

	int qType     = 0; // 0 ... Steinbach, 1 ... Sauter-Schwab
	int orderNear = 4; // 5 for Steinbach, 3-4 for Sauter-Schwab
	int orderFar  = 4; // 0 for none, else 3-4
};
#endif

namespace BEM4I {

bool isLinked()
{
#ifdef HAVE_BEM
	return true;
#else
	return false;
#endif
}

void deleteData(BEMData* &bem)
{
#ifdef HAVE_BEM
	bem4i::deleteBem4iData(bem->data);
	delete bem;
	bem = NULL;
#endif
}

void getLaplace(
		BEMData* &bem, double *K,
		esint nNodes, const double *nodes,
		esint nElements, const esint *elements,
		double conductivity)
{
#ifdef HAVE_BEM

	for (esint n = 0; n < nNodes; n++) {
		printf("%3.2f, %3.2f, %3.2f\n", nodes[3 * n + 0], nodes[3 * n + 1], nodes[3 * n + 2]);
	}
	for (esint n = 0; n < nElements; n++) {
		printf("%d %d %d\n", elements[3 * n + 0], elements[3 * n + 1], elements[3 * n + 2]);
	}
	printf("cond: %5.4f\n", conductivity);

	bem = new BEMData();
	bool nonsymmetric = false, verbose = false;
	bem4i::getLaplaceSteklovPoincare<esint, double>(
			K,
			NULL, NULL, // do not return the regularization matrix
			nNodes, nodes,
			nElements, elements,
			conductivity,
			bem->qType, bem->orderNear, bem->orderFar,
			bem->data,
			nonsymmetric, verbose);
#endif
}

void evaluateLaplace(
		BEMData* &bem, double *results,
		esint nNodes, const double *nodes,
		esint nElements, const esint *elements,
		esint nPoints, const double *points,
		double conductivity, double *dirichlet)
{
#ifdef HAVE_BEM

	bool verbose = false;
	bem4i::evaluateLaplaceRepresentationFormula<esint, double>(
			results,
			nNodes, nodes,
			nElements, elements,
			nPoints, points,
			conductivity, dirichlet,
			bem->orderFar,
			bem->data,
			verbose);
#endif
}

}
}
