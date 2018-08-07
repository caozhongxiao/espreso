
#ifndef SRC_INPUT_WORKBENCH_PARSER_CM_H_
#define SRC_INPUT_WORKBENCH_PARSER_CM_H_

#include "../parser/parser.h"

#include <map>

namespace espreso {

struct ESel;
struct PlainWorkbenchData;
struct MeshERegion;

struct NSel;
struct MeshNRegion;

struct CM: public WorkbenchParser {
	static size_t size;
	static const char* upper;
	static const char* lower;

	enum class Entity: int {
		VOLUME,
		AREA,
		LINE,
		KP,
		ELEMENTS,
		NODES
	};

	char name[MAX_NAME_SIZE];
	Entity entity;

	CM();
	CM& parse(const char* begin);

	bool addRegion(
			const PlainWorkbenchData &mesh,
			const std::vector<ESel> &esel, std::map<std::string, std::vector<eslocal> > &eregions,
			const std::vector<NSel> &nsel, std::map<std::string, std::vector<eslocal> > &nregions);

protected:
	bool addElementRegion(const PlainWorkbenchData &mesh, const std::vector<ESel> &esel, std::map<std::string, std::vector<eslocal> > &eregions);
	bool addNodeRegion(const std::vector<NSel> &nsel, std::map<std::string, std::vector<eslocal> > &nregions);
};

}



#endif /* SRC_INPUT_WORKBENCH_PARSER_CM_H_ */
