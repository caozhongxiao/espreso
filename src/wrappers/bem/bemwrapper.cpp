
#include "bemwrapper.h"

#include <cstddef>

#ifdef HAVE_BEM
#include "heatdtn.h"
#endif

namespace espreso {
namespace BEM4I {

struct BEMData
{
	bem4i::bem4idata<esint, double> * data;

	int qType     = 0; // 0 ... Steinbach, 1 ... Sauter-Schwab
	int orderNear = 4; // 5 for Steinbach, 3-4 for Sauter-Schwab
	int orderFar  = 4; // 0 for none, else 3-4
};

bool isLinked()
{
#ifdef HAVE_BEM
	return true;
#elif
	return false;
#endif
}

void deleteData(BEMData* &bem)
{
	bem4i::deleteBem4iData(bem->data);
}

void getLaplace(
		BEMData* &bem, double *K,
		esint nNodes, const double *nodes,
		esint nElements, const esint *elements,
		double conductivity)
{
#ifdef HAVE_BEM

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
