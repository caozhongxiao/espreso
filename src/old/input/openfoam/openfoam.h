
#ifndef INPUT_OPENFOAM_OPENFOAM_H_
#define INPUT_OPENFOAM_OPENFOAM_H_

#include "../loader.h"
#include "foam/boundary.h"
#include "foam/foamfile.h"
#include "foam/face.h"
#include "foam/dictionary.h"
#include "foam/elementbuilder.h"
#include "foam/zones.h"


namespace espreso {

struct InputConfiguration;

namespace input {

class OpenFOAM: public OldLoader {

public:
	static void load(const InputConfiguration &configuration, OldMesh &mesh, int rank, int size);
	bool faceBased() const { return true; }

protected:
	OpenFOAM(const InputConfiguration &configuration, OldMesh &mesh, int rank, int size);

	void points(Coordinates &coordinates);
	void elements(std::vector<size_t> &bodies, std::vector<OldElement*> &elements, std::vector<OldElement*> &faces, std::vector<OldElement*> &edges);
	void materials(std::vector<MaterialConfiguration*> &materials) {};
	void regions(
				std::vector<Evaluator*> &evaluators,
				std::vector<Region*> &regions,
				std::vector<OldElement*> &elements,
				std::vector<OldElement*> &faces,
				std::vector<OldElement*> &edges,
				std::vector<OldElement*> &nodes);
	void neighbours(std::vector<OldElement*> &nodes, std::vector<int> &neighbours, const std::vector<OldElement*> &faces, const std::vector<OldElement*> &edges);
	bool partitiate(const std::vector<OldElement*> &nodes, std::vector<eslocal> &partsPtrs, std::vector<std::vector<OldElement*> > &fixPoints, std::vector<OldElement*> &corners);

private:

	ParseError* computePolyMeshPath(int rank, int size);
	void solveParseError(ParseError *error);

	const InputConfiguration &_configuration;

	/** @brief Project path. */
	std::string _projectPath;

	/** @brief Path to PolyMesh, it contains also rank number for divided cases. */
	std::string _polyMeshPath;

	/** @brief Assigned rank, 0 for non MPI runs.*/
	int _rank;

	/** @brief Number of processes, 1 for non MPI runs*/
	int _size;
};

}
}

#endif /* INPUT_OPENFOAM_OPENFOAM_H_ */
