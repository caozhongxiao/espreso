
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

struct RegionMapBase;
class Mesh;
struct Element;
struct NodeData;
struct ProcessInterval;
struct BoundaryRegionStore;
template <typename TEBoundaries, typename TEData> class serializededata;
struct RBFTargetConfiguration;
struct RBFTargetTransformationConfiguration;
class Point;


class MeshPreprocessing {

	friend class Mesh;
public:

	void linkNodesAndElements();
	void exchangeHalo();

	void computeElementsNeighbors();
	void computeDecomposedDual(std::vector<eslocal> &dualDist, std::vector<eslocal> &dualData);

	void reclusterize(std::vector<eslocal> &dualDist, std::vector<eslocal> &dualData);
	void partitiate(eslocal parts, std::vector<eslocal> &dualDist, std::vector<eslocal> &dualData);

	void arrangeNodes();
	void arrangeElements();
	void arrangeRegions();

	void computeSharedFaceNodes();
	void computeCornerNodes();
	void computeFixPoints();
	void computeFixPointsOnSurface();
	void computeDomainsSurface();
	void triangularizeDomainSurface();

	void computeBodiesSurface();
	void computeSurfaceLocations();
	void searchContactInterfaces();

	void computeBoundaryElementsFromNodes(BoundaryRegionStore *bregion, int elementDimension);
	void computeRegionsIntersection(RegionMapBase &map);

	void morphRBF(const std::string &name, const RBFTargetConfiguration &configuration, int dimension);

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
	void exchangeElements(const std::vector<eslocal> &partition);
	void permuteElements(const std::vector<eslocal> &permutation, const std::vector<size_t> &distribution);
	void arrangeElementsPermutation(std::vector<eslocal> &permutation);
	void computeBoundaryNodes(std::vector<eslocal> &externalBoundary, std::vector<eslocal> &internalBoundary);
	void fillRegionMask();
	void computeRegionArea(BoundaryRegionStore *store);

	void addFixPoints(const serializededata<eslocal, eslocal>* elements, eslocal begin, eslocal end, const serializededata<eslocal, Element*>* epointers, std::vector<eslocal> &fixPoints);

	void synchronizeRegionNodes(const std::string &name, serializededata<eslocal, eslocal>* &rnodes, std::vector<ProcessInterval> &nintervals);
	void computeIntervalOffsets(std::vector<ProcessInterval> &intervals, eslocal &uniqueOffset, eslocal &uniqueSize, eslocal &uniqueTotalSize);

	void start(const std::string &message);
	void skip(const std::string &message);
	void finish(const std::string &message);

	void processMorpher(const RBFTargetTransformationConfiguration &target, int dimension,
		std::vector<Point> &sPoints, eslocal startPoint, std::vector<double> &sDisplacement);
	eslocal prepareMatrixM(std::vector<Point> &rPoints,
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
