
#include "cmblock.h"

#include "esinfo/mpiinfo.h"
#include "esinfo/envinfo.h"

#include "basis/containers/tarray.h"
#include "basis/utilities/parser.h"
#include "basis/logging/logging.h"

#include <algorithm>

using namespace espreso;

size_t CMBlock::size = 7;
const char* CMBlock::upper = "CMBLOCK";
const char* CMBlock::lower = "cmblock";

CMBlock::CMBlock()
: entity(Entity::NODE), NUMITEMS(-1),
  lineSize(-1), lineEndSize(-1),
  valueSize(-1), valueLength(-1)
{
	memset(name, '\0', MAX_NAME_SIZE);
}

CMBlock& CMBlock::parse(const char* begin)
{
	auto error = [&] (std::string &line) {
		ESINFO(ERROR) << "Workbench parse error: unknown format of CMBLOCK: " << line;
	};

	std::string commandLine = Parser::getLine(begin);
	std::string descriptionLine = Parser::getLine(begin + commandLine.size());

	lineEndSize = (*(commandLine.end() - 2) == '\r') ? 2 : 1;

	std::vector<std::string> command = Parser::split(commandLine, ",", false);
	std::vector<std::string> description = Parser::split(descriptionLine.substr(1, descriptionLine.size() - 2 - lineEndSize), "i");
	if (command.size() != 4) {
		error(commandLine);
	}
	if (description.size() != 2) {
		error(descriptionLine);
	}

	command[1] = Parser::strip(command[1]);
	memcpy(name, command[1].data(), command[1].size() < MAX_NAME_SIZE ? command[1].size() : MAX_NAME_SIZE);
	if (StringCompare::caseInsensitiveEq(command[2], "NODE")) {
		entity = Entity::NODE;
	} else {
		if (StringCompare::caseInsensitiveEq(command[2], "ELEMENT") || StringCompare::caseInsensitiveEq(command[2], "ELEM") || StringCompare::caseInsensitiveEq(command[2], "ELEMENTS")) {
			entity = Entity::ELEMENT;
		} else {
			error(commandLine);
		}
	}
	NUMITEMS = std::stol(command[3]);

	valueSize = std::stoi(description[0]);
	valueLength = std::stoi(description[1]);

	lineSize = valueSize * valueLength + lineEndSize;

	esint lastLineSize = 0;
	if ((NUMITEMS % valueSize)) {
		lastLineSize = (NUMITEMS % valueSize) * valueLength + lineEndSize;
	}

	WorkbenchParser::fillIndices(begin,
			begin + commandLine.size() + descriptionLine.size(),
			begin + commandLine.size() + descriptionLine.size() + (NUMITEMS / valueSize) * lineSize + lastLineSize);
	return *this;
}

bool CMBlock::readData(std::vector<esint> &indices)
{
	size_t threads = info::env::OMP_NUM_THREADS;

	const char *first = getFirst(), *last = getLast();
	esint size = (last - first) / lineSize;

	std::vector<esint> tdistribution = tarray<esint>::distribute(threads, size);
	std::vector<std::vector<esint> > tindices(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<char> value(valueLength + 1);
		for (auto data = first + lineSize * tdistribution[t]; data < first + lineSize * tdistribution[t + 1];) {
			for (esint n = 0; n < valueSize; ++n) {
				memcpy(value.data(), data, valueLength);
				data += valueLength;
				if (atol(value.data()) > 0) {
					tindices[t].push_back(atol(value.data()) - 1);
				}
			}
			data += lineEndSize;
		}
		if (lRank == info::mpi::MPIrank && t == threads - 1) {
			auto data = first + lineSize * tdistribution[t + 1];
			for (esint n = 0; n < NUMITEMS % valueSize; ++n) {
				memcpy(value.data(), data, valueLength);
				data += valueLength;
				if (atol(value.data()) > 0) {
					tindices[t].push_back(atol(value.data()) - 1);
				}
			}
		}
	}

	for (size_t t = 0; t < threads; t++) {
		indices.insert(indices.end(), tindices[t].begin(), tindices[t].end());
	}

	std::sort(indices.begin(), indices.end());

	return true;
}



