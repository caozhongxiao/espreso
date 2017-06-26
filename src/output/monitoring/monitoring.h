
#ifndef SRC_OUTPUT_MONITORING_MONITORING_H_
#define SRC_OUTPUT_MONITORING_MONITORING_H_

#include <fstream>

#include "../store.h"

namespace espreso {

class Mesh;
class Region;
enum class Property;
enum StatisticalData: int;

namespace output {

struct Monitor {
	size_t printSize;
	Region* region;
	Property property;
	StatisticalData statistics;
};

class Monitoring: public Store {

public:
	Monitoring(const OutputConfiguration &output, const Mesh *mesh);

	void storeSettings(const Step &step) {};
	void storeSettings(size_t steps) {};
	void storeSettings(const std::vector<size_t> &steps) {};

	void storeFETIData(const Step &step, const Instance &instance) {};

	void storeSolution(const Step &step, const std::vector<Solution*> &solution);
	void storeSolution(const Step &step, const std::vector<Solution*> &solution, const std::vector<std::pair<ElementType, Property> > &properties)
	{
		storeSolution(step, solution);
	}
	void finalize();

	static char delimiter;

protected:
	const Mesh *_mesh;
	std::ofstream _os;


	std::vector<Monitor> _monitors;
};

}
}


#endif /* SRC_OUTPUT_MONITORING_MONITORING_H_ */
