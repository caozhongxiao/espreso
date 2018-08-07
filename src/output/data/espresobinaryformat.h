
#ifndef SRC_OUTPUT_DATA_ESPRESOBINARYFORMAT_H_
#define SRC_OUTPUT_DATA_ESPRESOBINARYFORMAT_H_

namespace espreso {

class Mesh;
struct ECFRoot;

class ESPRESOBinaryFormat {

public:
	static void store(const Mesh &mesh, const ECFRoot &configuration);

};


}



#endif /* SRC_OUTPUT_DATA_ESPRESOBINARYFORMAT_H_ */
