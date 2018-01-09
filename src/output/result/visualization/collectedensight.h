
#ifndef SRC_OUTPUT_RESULT_VISUALIZATION_COLLECTEDENSIGHT_H_
#define SRC_OUTPUT_RESULT_VISUALIZATION_COLLECTEDENSIGHT_H_

#include <string>
#include <sstream>

#include "collectedvisualization.h"
#include "ensightwriter.h"

namespace espreso {

struct Step;
class Mesh;

struct CollectedEnSight: public CollectedVisualization {
	CollectedEnSight(const std::string &name, const Mesh &mesh);
	~CollectedEnSight();

	void updateMesh();
	void updateSolution(const Step &step);

protected:
	std::string codetotype(int code);
	void storecasefile();

	void storeDecomposition();

	std::string _path;
	std::string _name;

	std::stringstream _caseheader;
	std::stringstream _casegeometry;
	std::stringstream _casevariables;
	std::stringstream _casetime;

	int _variableCounter;

	const EnsightBinaryWriter _writer;
};

struct CollectedEnSightWithDecomposition: public virtual CollectedEnSight {
	CollectedEnSightWithDecomposition(const std::string &name, const Mesh &mesh): CollectedEnSight(name, mesh) {}

	void updateMesh()
	{
		CollectedEnSight::updateMesh();
		CollectedEnSight::storeDecomposition();
	}
};

}


#endif /* SRC_OUTPUT_RESULT_VISUALIZATION_COLLECTEDENSIGHT_H_ */
