
#ifndef SRC_OUTPUT_SOLUTION_VISUALIZATION_ENSIGHT_H_
#define SRC_OUTPUT_SOLUTION_VISUALIZATION_ENSIGHT_H_

#include "collectedvisualization.h"

#include <string>
#include <sstream>

namespace espreso {

struct Step;
class Mesh;

struct EnSight: public CollectedVisualization {
	EnSight(const std::string &name, const Mesh &mesh);
	~EnSight();

	void storeGeometry();
	void storeFETIData();
	void storeVariables(const Step &step);

protected:
	std::string codetotype(int code);
	void storecasefile();

	std::string _path;
	std::string _name;
	const Mesh &_mesh;

	std::stringstream _caseheader;
	std::stringstream _casegeometry;
	std::stringstream _casevariables;
	std::stringstream _casetime;

	int _variableCounter;
};

}



#endif /* SRC_OUTPUT_SOLUTION_VISUALIZATION_ENSIGHT_H_ */
