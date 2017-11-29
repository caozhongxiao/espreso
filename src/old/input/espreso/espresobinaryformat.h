
#ifndef INPUT_ESPRESO_ESPRESOBINARYFORMAT_H_
#define INPUT_ESPRESO_ESPRESOBINARYFORMAT_H_

#include "../loader.h"

namespace espreso {

struct InputConfiguration;

namespace input {

class ESPRESOBinaryFormat: public OldLoader {

public:
	static void load(const InputConfiguration &configuration, OldMesh &mesh, int rank, int size);

protected:
	ESPRESOBinaryFormat(const InputConfiguration &configuration, OldMesh &mesh, int rank, int size)
	: OldLoader(mesh), _configuration(configuration), _rank(rank), _size(size) { };

	void points(Coordinates &coordinates);
	void elements(std::vector<size_t> &bodies, std::vector<OldElement*> &elements, std::vector<OldElement*> &faces, std::vector<OldElement*> &edges);
	void materials(std::vector<MaterialConfiguration*> &materials);
	void regions(
			std::vector<OldEvaluator*> &evaluators,
			std::vector<Region*> &regions,
			std::vector<OldElement*> &elements,
			std::vector<OldElement*> &faces,
			std::vector<OldElement*> &edges,
			std::vector<OldElement*> &nodes);
	void neighbours(std::vector<OldElement*> &nodes, std::vector<int> &neighbours, const std::vector<OldElement*> &faces, const std::vector<OldElement*> &edges);
	bool partitiate(const std::vector<OldElement*> &nodes, std::vector<eslocal> &partsPtrs, std::vector<std::vector<OldElement*> > &fixPoints, std::vector<OldElement*> &corners);

private:
	const InputConfiguration &_configuration;
	int _rank;
	int _size;
};

}
}


#endif /* INPUT_ESPRESO_ESPRESOBINARYFORMAT_H_ */
