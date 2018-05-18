
#include "cm.h"

#include "esel.h"
#include "nsel.h"

#include "../../../basis/containers/tarray.h"
#include "../../../basis/utilities/parser.h"
#include "../../../basis/logging/logging.h"
#include "../../../config/ecf/environment.h"

#include <cstring>
#include <functional>
#include <algorithm>
#include "../../input.h"

using namespace espreso;

size_t CM::size = 3;
const char* CM::upper = "CM,";
const char* CM::lower = "cm,";

CM::CM()
: entity(Entity::ELEMENTS)
{
	memset(name, '\0', MAX_NAME_SIZE);
}

CM& CM::parse(const char* begin)
{
	std::string commandLine = Parser::getLine(begin);
	std::vector<std::string> command = Parser::split(Parser::strip(commandLine), ",", false);

	switch (command.size()) {
	case 3:
		if (StringCompare::caseInsensitiveEq("VOLU", command[2])) {
			entity = Entity::VOLUME;
			ESINFO(ERROR) << "ESPRESO Workbench parser error: implement cm, Entity='" << command[2] << "'";
		}
		if (StringCompare::caseInsensitiveEq("AREA", command[2])) {
			entity = Entity::AREA;
			ESINFO(ERROR) << "ESPRESO Workbench parser error: implement cm, Entity='" << command[2] << "'";
		}
		if (StringCompare::caseInsensitiveEq("LINE", command[2])) {
			entity = Entity::LINE;
			ESINFO(ERROR) << "ESPRESO Workbench parser error: implement cm, Entity='" << command[2] << "'";
		}
		if (StringCompare::caseInsensitiveEq("KP", command[2])) {
			entity = Entity::KP;
			ESINFO(ERROR) << "ESPRESO Workbench parser error: implement cm, Entity='" << command[2] << "'";
		}
		if (StringCompare::caseInsensitiveEq("ELEM", command[2])) {
			entity = Entity::ELEMENTS;
		}
		if (StringCompare::caseInsensitiveEq("NODE", command[2])) {
			entity = Entity::NODES;
			ESINFO(ERROR) << "ESPRESO Workbench parser error: implement cm, Entity='" << command[2] << "'";
		}
	case 2:
		memcpy(name, command[1].data(), command[1].size() < MAX_NAME_SIZE ? command[1].size() : MAX_NAME_SIZE);
		break;
	default:
		ESINFO(ERROR) << "ESPRESO Workbench parser error: unknown format of '" << commandLine << "'";
	}

	WorkbenchParser::fillIndices(begin, begin);
	return *this;
}

bool CM::addRegion(
		const std::vector<ESel> &esel, const std::vector<EData> &elements, std::vector<MeshERegion> &eregions,
		const std::vector<NSel> &nsel, std::vector<MeshNRegion> &nregions)
{
	switch (entity) {
	case Entity::ELEMENTS:
		return addElementRegion(esel, elements, eregions);
	case Entity::NODES:
		return addNodeRegion(nsel, nregions);
	default:
		return false;
	}
}

bool CM::addElementRegion(const std::vector<ESel> &esel, const std::vector<EData> &elements, std::vector<MeshERegion> &eregions)
{
	// clear relevant set means all elements
	std::vector<ESel> relevant;
	for (size_t i = 0; i < esel.size() && esel[i].first < first; i++) {
		switch (esel[i].type) {
		case ESel::Type::S:
		case ESel::Type::NONE:
			relevant.clear();
			relevant.push_back(esel[i]);
			break;
		case ESel::Type::A:
			if (relevant.size() == 1 && relevant[0].type == ESel::Type::NONE) {
				relevant.clear();
			}
			relevant.push_back(esel[i]);
			break;
		case ESel::Type::ALL:
			relevant.clear();
			break;
		default:
			ESINFO(ERROR) << "ESPRESO Workbench parser error: not implemented type of esel command '" << esel[i].command() << "'";
		}
	}

	if (relevant.size() == 0) {
		return true;
	}

	size_t threads = environment->OMP_NUM_THREADS;
	std::vector<size_t> edistribution = tarray<eslocal>::distribute(threads, elements.size());
	std::vector<std::vector<eslocal> > eid(threads);

	auto checkElements = [&] (std::function<void(size_t, eslocal)> chck) {
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (size_t e = edistribution[t]; e < edistribution[t + 1]; ++e) {
				chck(t, e);
			}
		}
	};

	eregions.push_back(MeshERegion());
	eregions.back().name = name;
	for (size_t i = 0; i < relevant.size(); i++) {
		switch (relevant[i].item) {

		case ESel::Item::TYPE:
			if (relevant[i].VMIN == relevant[i].VMAX) {
				checkElements([&] (size_t t, eslocal e) {
					if (elements[e].etype == relevant[i].VMIN) {
						eid[t].push_back(elements[e].id);
					}
				});
				break;
			}
			if (relevant[i].VINC == 1) {
				checkElements([&] (size_t t, eslocal e) {
					if (relevant[i].VMIN <= elements[e].etype && elements[e].etype <= relevant[i].VMAX) {
						eid[t].push_back(elements[e].id);
					}
				});
				break;
			}
			checkElements([&] (size_t t, eslocal e) {
				eslocal value = relevant[i].VMIN;
				if (relevant[i].VMIN <= elements[e].etype && elements[e].etype <= relevant[i].VMAX && (elements[e].etype - relevant[i].VMIN) % relevant[i].VINC == 0) {
					eid[t].push_back(elements[e].id);
				}
			});
			break;

		case ESel::Item::MAT:
			if (relevant[i].VMIN == relevant[i].VMAX) {
				checkElements([&] (size_t t, eslocal e) {
					if (elements[e].material == relevant[i].VMIN) {
						eid[t].push_back(elements[e].id);
					}
				});
				break;
			}
			if (relevant[i].VINC == 1) {
				checkElements([&] (size_t t, eslocal e) {
					if (relevant[i].VMIN <= elements[e].material && elements[e].material <= relevant[i].VMAX) {
						eid[t].push_back(elements[e].id);
					}
				});
				break;
			}
			checkElements([&] (size_t t, eslocal e) {
				eslocal value = relevant[i].VMIN;
				if (relevant[i].VMIN <= elements[e].material && elements[e].material <= relevant[i].VMAX && (elements[e].material - relevant[i].VMIN) % relevant[i].VINC == 0) {
					eid[t].push_back(elements[e].id);
				}
			});
			break;

		default:
			ESINFO(ERROR) << "ESPRESO Workbench parser error: not implemented type of esel command '" << esel[i].command() << "'";
		}
	}

	for (size_t t = 0; t < threads; t++) {
		eregions.back().elements.insert(eregions.back().elements.end(), eid[t].begin(), eid[t].end());
	}
	std::sort(eregions.back().elements.begin(), eregions.back().elements.end());

	return true;
}

bool CM::addNodeRegion(const std::vector<NSel> &nsel, std::vector<MeshNRegion> &nregions)
{
	return false;
}
