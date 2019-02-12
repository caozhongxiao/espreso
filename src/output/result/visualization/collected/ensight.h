
#ifndef SRC_OUTPUT_RESULT_VISUALIZATION_COLLECTED_ENSIGHT_H_
#define SRC_OUTPUT_RESULT_VISUALIZATION_COLLECTED_ENSIGHT_H_

#include <string>
#include <sstream>

#include "collectedvisualization.h"
#include "output/result/visualization/ensightwriter.h"

namespace espreso {

class Mesh;

struct EnSight: public CollectedVisualization {
	EnSight(const std::string &name, const Mesh &mesh);
	~EnSight();

	void updateMesh();
	void updateSolution();

protected:
	std::string codetotype(int code);
	void storecasefile();
	void setvariables();

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

struct EnSightWithDecomposition: public virtual EnSight {
	EnSightWithDecomposition(const std::string &name, const Mesh &mesh): EnSight(name, mesh) {}

	void updateMesh()
	{
		EnSight::updateMesh();
		EnSight::storeDecomposition();
	}
};

}


#endif /* SRC_OUTPUT_RESULT_VISUALIZATION_COLLECTED_ENSIGHT_H_ */
