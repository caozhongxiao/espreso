//
//#ifndef SRC_MESH_TRANSFORMATION_TRANSFORMATIONS_H_
//#define SRC_MESH_TRANSFORMATION_TRANSFORMATIONS_H_
//
//#include "tflags.h"
//
//#include <cstddef>
//#include <vector>
//
//namespace espreso {
//
//class Mesh;
//class Element;
//class ElementStore;
//class BoundaryStore;
//class EInterval;
//template <typename TEBoundaries, typename TEData> class serializededata;
//
//#define TVERBOSITY PROGRESS1
//
//struct Transformation {
//
//	// static void addLinkFromTo(Mesh &mesh, TFlags::ELEVEL from, TFlags::ELEVEL to);
//
//	// static void exchangeHaloElements(Mesh &mesh);
//
////	static void computeDual(Mesh &mesh);
////	static void computeDecomposedDual(Mesh &mesh, TFlags::SEPARATE separate);
//
//	static void computeProcessBoundaries(Mesh &mesh);
//	static void computeDomainsBoundaries(Mesh &mesh);
//
////	static void computeElementCenters(Mesh &mesh);
////	static void computeDomainsCenters(Mesh &mesh);
//
////	static void reclusterize(Mesh &mesh);
////	static void partitiate(Mesh &mesh, esglobal parts, TFlags::SEPARATE separate);
//
////	static void reindexNodes(Mesh &mesh);
////	static void arrangeNodes(Mesh &mesh);
//
////	static void assignDomainsToNodes(Mesh &mesh);
//private:
//	static size_t level;
//
////	static void exchangeElements(Mesh &mesh, const std::vector<esglobal> &partition);
////	static void permuteElements(Mesh &mesh, const std::vector<eslocal> &permutation, const std::vector<size_t> &distribution);
////	static void projectNodesToDomains(Mesh &mesh);
//	static void computeIntervals(std::vector<EInterval> &intervals, const serializededata<eslocal, eslocal> &compdata, const std::vector<size_t> &distribution, const std::vector<eslocal> &permutation);
//
//	template <typename Tdual>
//	static void computeBoundaries(Mesh &mesh,
//			serializededata<eslocal, Tdual>       *elementDual,
//			esglobal                               dualOffset,
//			const std::vector<esglobal>            &IDBoundaries,
//			std::vector<std::vector<eslocal> >     *elementData,
//			std::vector<std::vector<eslocal> >     *faceDistribution,
//			std::vector<std::vector<eslocal> >     *faceData,
//			std::vector<std::vector<Element*> > *faceCodes,
//			std::vector<std::vector<int> >         *faceNeighbors);
//
//	static void distributeElementsToIntervals(Mesh &mesh,
//			BoundaryStore*                         &boundaries,
//			const std::vector<esglobal>            &IDBoundaries,
//			std::vector<std::vector<eslocal> >     &elementData);
//
//	static void distributeFacesToIntervals(Mesh &mesh,
//			BoundaryStore*                         &boundaries,
//			std::vector<std::vector<eslocal> >     &faceDistribution,
//			std::vector<std::vector<eslocal> >     &faceData,
//			std::vector<std::vector<Element*> > &faceCodes,
//			std::vector<std::vector<int> >         &faceNeighbors);
//
//	static void distributeNodesToIntervals(Mesh &mesh, BoundaryStore* &boundaries);
//};
//
//}
//
//#endif /* SRC_MESH_TRANSFORMATION_TRANSFORMATIONS_H_ */
