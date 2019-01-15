
#ifndef SRC_OUTPUT_RESULT_VISUALIZATION_COLLECTED_STL_H_
#define SRC_OUTPUT_RESULT_VISUALIZATION_COLLECTED_STL_H_

#include <string>
#include <sstream>

#include "collectedvisualization.h"
#include "output/result/visualization/stlwritter.h"

namespace espreso {

class Mesh;

struct STL: public CollectedVisualization {
	STL(const std::string &name, const Mesh &mesh);
	~STL();

	void updateMesh();
	void updateSolution();

protected:
	std::string _path;
	std::string _name;

	const STLBinaryWriter _writer;
};

}


#endif /* SRC_OUTPUT_RESULT_VISUALIZATION_COLLECTED_STL_H_ */
