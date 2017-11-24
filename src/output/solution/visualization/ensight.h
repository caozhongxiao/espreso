
#ifndef SRC_OUTPUT_SOLUTION_VISUALIZATION_ENSIGHT_H_
#define SRC_OUTPUT_SOLUTION_VISUALIZATION_ENSIGHT_H_

#include <string>

namespace espreso {

class Mesh;

struct EnSight {
	EnSight(const std::string &name, const Mesh &mesh);
	~EnSight();

	void storeGeometry();
	void storeVariables();

protected:
	std::string _path;
	std::string _name;
	const Mesh &_mesh;

	std::ofstream *_casefile;
};

}



#endif /* SRC_OUTPUT_SOLUTION_VISUALIZATION_ENSIGHT_H_ */
