
#ifndef SRC_MESH_TRANSFORMATION_TRANSFORMATIONS_H_
#define SRC_MESH_TRANSFORMATION_TRANSFORMATIONS_H_

#include "tflags.h"

#include <cstddef>
#include <vector>

namespace espreso {

class NewMesh;
class NewElement;
class ElementStore;
class BoundaryStore;
class EInterval;
template <typename TEBoundaries, typename TEData> class serializededata;

#define TVERBOSITY PROGRESS1

struct Transformation {

	static void addLinkFromTo(NewMesh &mesh, TFlags::ELEVEL from, TFlags::ELEVEL to);

	static void exchangeHaloElements(NewMesh &mesh);

	static void computeDual(NewMesh &mesh);
	static void computeDecomposedDual(NewMesh &mesh, TFlags::SEPARATE separate);

	static void computeProcessBoundaries(NewMesh &mesh);
	static void computeDomainsBoundaries(NewMesh &mesh);

	static void computeElementCenters(NewMesh &mesh);
	static void computeDomainsCenters(NewMesh &mesh);

	static void reclusterize(NewMesh &mesh);
	static void partitiate(NewMesh &mesh, esglobal parts, TFlags::SEPARATE separate);

	static void reindexNodes(NewMesh &mesh);
	static void arrangeNodes(NewMesh &mesh);

	static void assignDomainsToNodes(NewMesh &mesh);
private:
	static size_t level;

	static void exchangeElements(NewMesh &mesh, const std::vector<esglobal> &partition);
	static void permuteElements(NewMesh &mesh, const std::vector<eslocal> &permutation, const std::vector<size_t> &distribution);
	static void projectNodesToDomains(NewMesh &mesh);
	static void computeIntervals(std::vector<EInterval> &intervals, const serializededata<eslocal, eslocal> &compdata, const std::vector<size_t> &distribution, const std::vector<eslocal> &permutation);

	template <typename Tdual>
	static void computeBoundaries(NewMesh &mesh,
			serializededata<eslocal, Tdual>       *elementDual,
			esglobal                               dualOffset,
			const std::vector<esglobal>            &IDBoundaries,
			std::vector<std::vector<eslocal> >     *elementData,
			std::vector<std::vector<eslocal> >     *faceDistribution,
			std::vector<std::vector<eslocal> >     *faceData,
			std::vector<std::vector<NewElement*> > *faceCodes,
			std::vector<std::vector<int> >         *faceNeighbors);

	static void distributeElementsToIntervals(NewMesh &mesh,
			BoundaryStore*                         &boundaries,
			const std::vector<esglobal>            &IDBoundaries,
			std::vector<std::vector<eslocal> >     &elementData);

	static void distributeFacesToIntervals(NewMesh &mesh,
			BoundaryStore*                         &boundaries,
			std::vector<std::vector<eslocal> >     &faceDistribution,
			std::vector<std::vector<eslocal> >     &faceData,
			std::vector<std::vector<NewElement*> > &faceCodes,
			std::vector<std::vector<int> >         &faceNeighbors);

	static void distributeNodesToIntervals(NewMesh &mesh, BoundaryStore* &boundaries);
};

}



#endif /* SRC_MESH_TRANSFORMATION_TRANSFORMATIONS_H_ */
