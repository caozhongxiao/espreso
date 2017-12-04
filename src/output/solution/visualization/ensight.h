
#ifndef SRC_OUTPUT_SOLUTION_VISUALIZATION_ENSIGHT_H_
#define SRC_OUTPUT_SOLUTION_VISUALIZATION_ENSIGHT_H_

#include "collectedvisualization.h"

#include <string>

namespace espreso {

class Mesh;

struct EnSight: public CollectedVisualization {
	EnSight(const std::string &name, const Mesh &mesh);
	~EnSight();

	void storeGeometry();
	void storeFETIData();
	void storeVariables();

protected:
	std::string codetotype(int code);

	std::string _path;
	std::string _name;
	const Mesh &_mesh;

	std::ofstream *_casefile;
};

}



#endif /* SRC_OUTPUT_SOLUTION_VISUALIZATION_ENSIGHT_H_ */
