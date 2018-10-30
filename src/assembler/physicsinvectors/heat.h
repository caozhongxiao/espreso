
#ifndef SRC_ASSEMBLER_PHYSICSINVECTORS_HEAT_H_
#define SRC_ASSEMBLER_PHYSICSINVECTORS_HEAT_H_

#include "physicsinvectors.h"

namespace espreso {

enum class Property;
class MaterialBaseConfiguration;
class HeatTransferConfiguration;
class NodeData;

struct Heat: public PhysicsInVectors
{
	Heat(Mesh &mesh, Instance &instance, Step &step, const HeatTransferConfiguration &configuration);

	void initData();
	void updateData();

	void processElement(eslocal domain, eslocal eindex, DenseMatrix &Ke, DenseMatrix &fe);

	// node values
	serializededata<eslocal, double> *_coordinates, *_K;

	// element values
	serializededata<eslocal, double>* _area;

	NodeData *_temperature;
};

}



#endif /* SRC_ASSEMBLER_PHYSICSINVECTORS_HEAT_H_ */
