
#include "nblock.h"

#include "../../../basis/containers/point.h"
#include "../../../basis/containers/tarray.h"
#include "../../../basis/utilities/parser.h"
#include "../../../basis/utilities/utils.h"

#include "../../../config/ecf/environment.h"

using namespace espreso;

size_t NBlock::size = 6;
const char* NBlock::upper = "NBLOCK";
const char* NBlock::lower = "nblock";

NBlock::NBlock()
: NUMFIELD(-1), Solkey(0), NDMAX(-1), NDSEL(-1),
  lineSize(-1), lineEndSize(-1),
  indexSize(-1), indexLength(-1), valueSize(-1), valueLength(-1)
{

}

NBlock& NBlock::parse(const char* begin)
{
	auto error = [&] (std::string &line) {
		ESINFO(ERROR) << "Workbench parse error: unknown format of NBLOCK: " << line;
	};

	std::string commandLine = Parser::getLine(begin);
	std::string descriptionLine = Parser::getLine(begin + commandLine.size());

	lineEndSize = (*(commandLine.end() - 2) == '\r') ? 2 : 1;

	std::vector<std::string> command = Parser::split(commandLine, ",", false);
	std::vector<std::string> description = Parser::split(descriptionLine.substr(1, descriptionLine.size() - 2 - lineEndSize), ",");
	if (description.size() != 2) {
		error(descriptionLine);
	}
	std::vector<std::string> indices = Parser::split(description[0], "i");
	std::vector<std::string> values = Parser::split(description[1], "e");

	if (indices.size() != 2 || values.size() < 2) {
		error(descriptionLine);
	}

	switch (command.size()) {
	case 5:
		NDSEL = command[4].size() ? std::stol(command[4]) : -1;
	case 4:
		NDMAX = command[3].size() ? std::stol(command[3]) : -1;
	case 3:
		Solkey = command[2].size();
	case 2:
		NUMFIELD = std::stoi(command[1]);
		break;
	default:
		error(commandLine);
	}

	indexSize = std::stoi(indices[0]);
	indexLength = std::stoi(indices[1]);
	valueSize = std::stoi(values[0]);
	valueLength = std::stoi(values[1]);

	lineSize = indexSize * indexLength + 3 * valueLength + lineEndSize; // rotation can be omitted

	WorkbenchParser::fillIndices(begin, begin + commandLine.size() + descriptionLine.size());
	return *this;
}

bool NBlock::readData(std::vector<eslocal> &nIDs, std::vector<Point> &coordinates, double scaleFactor)
{
	if (Solkey) {
		if (NUMFIELD == 6) {
			return index_solid_line_x_y_z(nIDs, coordinates, scaleFactor);
		}
	} else {
		if (NUMFIELD == 3) {
			return index_x_y_z(nIDs, coordinates, scaleFactor);
		}
	}
	return false;
}

bool NBlock::index_x_y_z(std::vector<eslocal> &nIDs, std::vector<Point> &coordinates, double scaleFactor)
{
	size_t threads = environment->OMP_NUM_THREADS;

	const char *first = getFirst(), *last = getLast();
	size_t size = (last - first) / lineSize;

	size_t offset = coordinates.size();
	nIDs.resize(offset + size);
	coordinates.resize(offset + size);
	std::vector<size_t> tdistribution = tarray<size_t>::distribute(threads, size);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<char> value(valueLength + 1);
		std::vector<char> index(indexLength + 1);
		eslocal i = 0;
		double x, y, z;
		for (auto data = first + lineSize * tdistribution[t]; data < first + lineSize * tdistribution[t + 1]; ++i) {
			memcpy(index.data(), data, indexLength);
			data += indexLength;
			nIDs[offset + tdistribution[t] + i] = atol(index.data()) - 1;

			memcpy(value.data(), data, valueLength);
			data += valueLength;
			x = scaleFactor * atof(value.data());

			memcpy(value.data(), data, valueLength);
			data += valueLength;
			y = scaleFactor * atof(value.data());

			memcpy(value.data(), data, valueLength);
			data += valueLength;
			z = scaleFactor * atof(value.data());

			data += lineEndSize;

			coordinates[offset + tdistribution[t] + i] = Point(x, y, z);
		}
	}
	return true;
}

bool NBlock::index_solid_line_x_y_z(std::vector<eslocal> &nIDs, std::vector<Point> &coordinates, double scaleFactor)
{
	size_t threads = environment->OMP_NUM_THREADS;

	const char *first = getFirst(), *last = getLast();
	eslocal size = (last - first) / lineSize;

	size_t offset = coordinates.size();
	nIDs.resize(offset + size);
	coordinates.resize(offset + size);
	std::vector<size_t> tdistribution = tarray<eslocal>::distribute(threads, size);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<char> value(valueLength + 1);
		std::vector<char> index(indexLength + 1);
		eslocal i = 0;
		double x, y, z;
		for (auto data = first + lineSize * tdistribution[t]; data < first + lineSize * tdistribution[t + 1]; ++i) {
			memcpy(index.data(), data, indexLength);
			data += indexLength;
			nIDs[offset + tdistribution[t] + i] = atol(index.data()) - 1;

			data += (indexSize - 1)* indexLength; // skip index, solid, line

			memcpy(value.data(), data, valueLength);
			data += valueLength;
			x = scaleFactor * atof(value.data());

			memcpy(value.data(), data, valueLength);
			data += valueLength;
			y = scaleFactor * atof(value.data());

			memcpy(value.data(), data, valueLength);
			data += valueLength;
			z = scaleFactor * atof(value.data());

			data += lineEndSize;

			coordinates[offset + tdistribution[t] + i] = Point(x, y, z);
		}
	}
	return true;
}
