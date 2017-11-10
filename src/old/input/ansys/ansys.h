
#ifndef INPUT_ANSYS_ANSYS_H_
#define INPUT_ANSYS_ANSYS_H_

#include <string>
#include <sstream>
#include <fstream>
#include <vector>

#include "../loader.h"
#include "utils.h"
#include "parser.h"

namespace espreso {

struct InputConfiguration;

namespace input {

class AnsysWorkbench: public OldLoader {

public:
	static void load(const InputConfiguration &configuration, Mesh &mesh, int rank, int size);

protected:
	AnsysWorkbench(const InputConfiguration &configuration, Mesh &mesh, int rank, int size)
	: OldLoader(mesh), _workbench(configuration), _parser(mesh) { };

	void points(Coordinates &coordinates);
	void elements(std::vector<size_t> &bodies, std::vector<Element*> &elements, std::vector<Element*> &faces, std::vector<Element*> &edges);
	void materials(std::vector<MaterialConfiguration*> &materials);
	void regions(
			std::vector<Evaluator*> &evaluators,
			std::vector<Region*> &regions,
			std::vector<Element*> &elements,
			std::vector<Element*> &faces,
			std::vector<Element*> &edges,
			std::vector<Element*> &nodes);
	void neighbours(std::vector<Element*> &nodes, std::vector<int> &neighbours, const std::vector<Element*> &faces, const std::vector<Element*> &edges);
	bool partitiate(const std::vector<Element*> &nodes, std::vector<eslocal> &partsPtrs, std::vector<std::vector<Element*> > &fixPoints, std::vector<Element*> &corners);

	void open();
	void close();

private:
	const InputConfiguration &_workbench;
	WorkbenchParser _parser;
};

}
}





#endif /* INPUT_ANSYS_ANSYS_H_ */
