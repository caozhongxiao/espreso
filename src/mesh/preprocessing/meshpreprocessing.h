
#ifndef SRC_MESH_PREPROCESSING_MESHPREPROCESSING_H_
#define SRC_MESH_PREPROCESSING_MESHPREPROCESSING_H_

#include <cstddef>
#include <string>
#include <vector>

namespace espreso {

#define VERBOSITY(level) level == 0 ? PROGRESS2 : PROGRESS3

class Mesh;
struct Element;
struct ProcessInterval;
struct BoundaryRegionStore;
template <typename TEBoundaries, typename TEData> class serializededata;

class MeshPreprocessing {

	friend class Mesh;
public:

	void linkNodesAndElements();
	void exchangeHalo();

	void computeDual();
	void computeDecomposedDual(bool separateMaterials, bool separateRegions, bool separateEtype);
	void computeFullDual(const serializededata<eslocal, eslocal>* elements, eslocal begin, eslocal end, std::vector<eslocal> &dist, std::vector<eslocal> &data);

	void reclusterize();
	void partitiate(eslocal parts, bool separateMaterials, bool separateRegions, bool separateEtype);

	void arrangeNodes();
	void arrangeElements();
	void arrangeRegions();

	void computeSharedFaceNodes();
	void computeCornerNodes();
	void computeFixPoints();
	void computeFixPointsOnSurface();
	void computeDomainsSurface();
	void triangularizeDomainSurface();

protected:
	static size_t level;

	MeshPreprocessing(Mesh *mesh): _mesh(mesh) {}
	Mesh *_mesh;

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
};

}



#endif /* SRC_MESH_PREPROCESSING_MESHPREPROCESSING_H_ */
