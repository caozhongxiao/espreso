
#ifndef INPUT_ESPRESO_ESPRESOBINARYFORMAT_H_
#define INPUT_ESPRESO_ESPRESOBINARYFORMAT_H_

#include "../loader.h"

namespace espreso {

struct InputConfiguration;

namespace input {

class ESPRESOBinaryFormat: public OldLoader {

public:
	static void load(const InputConfiguration &configuration, Mesh &mesh, int rank, int size);

protected:
	ESPRESOBinaryFormat(const InputConfiguration &configuration, Mesh &mesh, int rank, int size)
	: OldLoader(mesh), _configuration(configuration), _rank(rank), _size(size) { };

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

private:
	const InputConfiguration &_configuration;
	int _rank;
	int _size;
};

}
}


#endif /* INPUT_ESPRESO_ESPRESOBINARYFORMAT_H_ */
