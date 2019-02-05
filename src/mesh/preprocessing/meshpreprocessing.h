
#ifndef SRC_MESH_PREPROCESSING_MESHPREPROCESSING_H_
#define SRC_MESH_PREPROCESSING_MESHPREPROCESSING_H_

#include <cstddef>
#include <string>
#include <vector>
#include <map>

namespace espreso {

#define VERBOSITY(level) level == 0 ? PROGRESS2 : PROGRESS3

class TimeEval;
class TimeEvent;

class Mesh;
struct Element;
struct NodeData;
struct ProcessInterval;
struct BoundaryRegionStore;
struct SurfaceStore;
template <typename TEBoundaries, typename TEData> class serializededata;
struct RBFTargetConfiguration;
struct RBFTargetTransformationConfiguration;
class Point;


class MeshPreprocessing {

	friend class Mesh;
public:

	void linkNodesAndElements();
	void exchangeHalo();
	void exchangeElements(const std::vector<esint> &partition);

	void computeElementsNeighbors();
	void computeElementsCenters();
	void computeDecomposedDual(std::vector<esint> &dualDist, std::vector<esint> &dualData);

	void reclusterize();
	void partitiate(esint parts, bool uniformDecomposition);

	void arrangeNodes();
	void arrangeElements();
	void arrangeRegions();

	void computeLocalIndices();
	void computeSharedFaceNodes();
	void computeCornerNodes();
	void computeFixPoints();
	void computeFixPointsOnSurface();
	void computeDomainsSurface();
	void triangularizeDomainSurface();

	void triangularizeSurface(SurfaceStore *surface);
	void triangularizeBoundary(BoundaryRegionStore *boundary);

	void computeRegionsSurface();
	void computeBodiesSurface();
	void computeSurfaceLocations();
	void searchContactInterfaces();

	void computeBoundaryElementsFromNodes(BoundaryRegionStore *bregion, int elementDimension);

	void morphRBF(const std::string &name, const RBFTargetConfiguration &configuration, int dimension);

	void startPreprocessing();
	void finishPreprocessing();

protected:
	static size_t level;

	MeshPreprocessing(Mesh *mesh);
	~MeshPreprocessing();
	Mesh *_mesh;

	NodeData* _morphing;

	TimeEval *_timeStatistics;
	std::map<std::string, TimeEvent*> _timeEvents;

private:
	void permuteElements(const std::vector<esint> &permutation, const std::vector<size_t> &distribution);
	void arrangeElementsPermutation(std::vector<esint> &permutation);
	void computeBoundaryNodes(std::vector<esint> &externalBoundary, std::vector<esint> &internalBoundary);
	void fillRegionMask();
	void computeRegionArea(BoundaryRegionStore *store);

	void addFixPoints(const serializededata<esint, esint>* elements, esint begin, esint end, const serializededata<esint, Element*>* epointers, std::vector<esint> &fixPoints);

	void synchronizeRegionNodes(const std::string &name, serializededata<esint, esint>* &rnodes, std::vector<ProcessInterval> &nintervals);
	void computeIntervalOffsets(std::vector<ProcessInterval> &intervals, esint &uniqueOffset, esint &uniqueSize, esint &uniqueTotalSize);

	void start(const std::string &message);
	void skip(const std::string &message);
	void finish(const std::string &message);

	void processMorpher(const RBFTargetTransformationConfiguration &target, int dimension,
		std::vector<Point> &sPoints, esint startPoint, std::vector<double> &sDisplacement);
	esint prepareMatrixM(std::vector<Point> &rPoints,
			std::vector<double> &rDisplacement,
			int dimension, const RBFTargetConfiguration &configuration,
			std::vector<double> &M_values,
			bool use_x = true,bool use_y = true,bool use_z = true);
	void readExternalFile(
			const RBFTargetConfiguration &configuration, int dimension,
			std::map<std::string, std::vector<Point>> &external_data);
};

}



#endif /* SRC_MESH_PREPROCESSING_MESHPREPROCESSING_H_ */
