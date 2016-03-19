
#ifndef INPUT_ANSYS_MATSOL_ANSYS_H_
#define INPUT_ANSYS_MATSOL_ANSYS_H_

#include <string>
#include <sstream>
#include <fstream>

#include "../../loader.h"
#include "../utils.h"

namespace espreso {
namespace input {

class AnsysMatsol: public Loader {

public:
	static void load(Mesh &mesh, const Options &options, int rank, int size)
	{
		AnsysMatsol matsol(mesh, options, rank, size);
		matsol.fill();
	}

protected:
	AnsysMatsol(Mesh &mesh, const Options &options, int rank, int size)
	: Loader(mesh), _path(options.path) { };

	void points(Coordinates &coordinates);
	void elements(std::vector<Element*> &elements);
	void boundaryConditions(Coordinates &coordinates);
	void clusterBoundaries(Boundaries &boundaries, std::vector<int> &neighbours);

private:
	static size_t getLinesCount(const std::string &file);
	struct TestEOL {
		bool operator()(char c) {
			return c == '\n';
		}
	};

	std::string _path;
};

}
}




#endif /* INPUT_ANSYS_MATSOL_ANSYS_H_ */
