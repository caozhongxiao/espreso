
#include "workbench.h"

#include "../../basis/containers/tarray.h"
#include "../../basis/logging/logging.h"
#include "../../basis/utilities/parser.h"
#include "../../basis/utilities/communication.h"
#include "../../basis/utilities/utils.h"

#include "../../config/ecf/ecf.h"

#include <iostream>
#include <string>
#include <regex>

#define MAX_LINE_SIZE 500 // upper bound on line size
#define MAX_LINE_STEP 2   // sometimes we need to red more lines to get full information

using namespace espreso;

void WorkbenchLoader::load(const ECFConfiguration &configuration, Mesh &mesh)
{
	WorkbenchLoader(configuration, mesh);
}

WorkbenchLoader::WorkbenchLoader(const ECFConfiguration &configuration, Mesh &mesh)
: _configuration(configuration), _mesh(mesh)
{
	readData();
	parseData();
}

int WorkbenchLoader::elementNodeCount(int etype)
{
	switch (etype) {
	case 186:
		return 20;
	default:
		ESINFO(ERROR) << "Workbench parser error: not implemented element type " << etype;
		return -1;
	}
}

void WorkbenchLoader::readData()
{
	MPI_File MPIfile;

	if (MPI_File_open(environment->MPICommunicator, _configuration.workbench.path.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &MPIfile)) {
		ESINFO(ERROR) << "MPI cannot create file '" << _configuration.workbench.path << "'";
	}

	MPI_Offset size;
	MPI_File_get_size(MPIfile, &size);

	MPI_Datatype fDataDistribution;
	std::vector<size_t> fdistribution = tarray<int>::distribute(environment->MPIsize, size);

	std::vector<MPI_Aint> displacement;
	std::vector<int> length;

	displacement.push_back(fdistribution[environment->MPIrank]);
	length.push_back(fdistribution[environment->MPIrank + 1] - fdistribution[environment->MPIrank]);

	MPI_Type_create_hindexed(1, length.data(), displacement.data(), MPI_BYTE, &fDataDistribution);
	MPI_Type_commit(&fDataDistribution);

	_data.resize(length.front() + MAX_LINE_STEP * MAX_LINE_SIZE);

	MPI_File_set_view(MPIfile, 0, MPI_BYTE, fDataDistribution, "native", MPI_INFO_NULL);
	MPI_File_read_all(MPIfile, _data.data(), _data.size() - MAX_LINE_STEP * MAX_LINE_SIZE, MPI_BYTE, MPI_STATUS_IGNORE);
	MPI_File_close(&MPIfile);

	_current = _data.data();
	if (environment->MPIsize > 1) { // align to line end, TODO: fix for tiny files
		char *exchangeStart = _current;
		if (environment->MPIrank) {
			while (*_current++ != '\n');
			exchangeStart = _current;
			for (int i = 1; i < MAX_LINE_STEP; i++) {
				while (*exchangeStart++ != '\n');
			}
		}

		std::vector<std::vector<char> > sBuffer(1, std::vector<char>(_data.data(), exchangeStart)), rBuffer;

		std::vector<int> neighs;
		if (environment->MPIrank) {
			neighs.push_back(environment->MPIrank - 1);
		}
		if (environment->MPIrank + 1 < environment->MPIsize) {
			neighs.push_back(environment->MPIrank + 1);
		}
		rBuffer.resize(neighs.size());
		Communication::receiveUpperUnknownSize(sBuffer, rBuffer, neighs);

		if (rBuffer.back().size() > MAX_LINE_STEP * MAX_LINE_SIZE) {
			ESINFO(ERROR) << "ESPRESO internal error: increase max line size in Ansys file.";
		}
		memcpy(_data.data() + _data.size() - MAX_LINE_STEP * MAX_LINE_SIZE, rBuffer.back().data(), rBuffer.back().size());
		char* firstLineEnd = rBuffer.back().data();
		if (rBuffer.back().size()) {
			while (*firstLineEnd++ != '\n');
		}
		_end = _data.data() + _data.size() - MAX_LINE_STEP * MAX_LINE_SIZE + (firstLineEnd - rBuffer.back().data());
	} else {
		_end = _data.data() + _data.size() - MAX_LINE_STEP * MAX_LINE_SIZE;
	}
	_begin = _current;
	_dataOffset = Communication::getDistribution<eslocal>(_end - _current);
}

std::string WorkbenchLoader::getLine(eslocal index)
{
	char *first = _begin + index - _dataOffset[environment->MPIrank];
	char *last = first;
	while (*last != '\n' && *last != '\r') { ++last; }
	return std::string(first, last);
}

void WorkbenchLoader::parseData()
{
	auto equal = [] (const std::string &param, const std::cmatch &match) {
		return param.size() == match.length(0) && std::equal(param.begin(), param.end(), match[0].first);
	};

	auto mindex = [&] () {
		return _dataOffset[environment->MPIrank] + _current - _begin;
	};

	auto fillranks = [&] (DataInterval &data) {
		data.sRank = std::lower_bound(_dataOffset.begin(), _dataOffset.end(), data.sIndex + 1) - _dataOffset.begin() - 1;
		data.eRank = std::lower_bound(_dataOffset.begin(), _dataOffset.end(), data.eIndex + 1) - _dataOffset.begin() - 1;
	};

	std::vector<eslocal> blockends, etypes;

	std::string coordinates = "nblock,";
	std::string blockend    = "\\s-1\\s";
	std::string etype       = "\\set,";
	std::string eblock      = "eblock,";
	std::string cmblock     = "CMBLOCK,";


	std::regex regex(
			coordinates + "|" + blockend + "|" + etype + "|" + eblock + "|" + cmblock,
			std::regex::optimize);

	std::cmatch match;

	while (std::regex_search(_current, match, regex) && match.position(0) < _end - _current) {
		_current += match.position(0);

		if (equal(coordinates, match)) {
			_coordinates.push_back(DataInterval());
			_coordinates.back().header = mindex();
		}
		if (equal(eblock, match)) {
			_elements.push_back(DataInterval());
			_elements.back().header = mindex();
		}
		if (equal(cmblock, match)) {
			_cmblocks.push_back(DataInterval());
			_cmblocks.back().header = mindex();
		}

		if (match[0].str().find("et") != -1) {
			etypes.push_back(mindex() + match[0].str().find("et"));
		}
		if (match[0].str().find("-1") != -1) {
			blockends.push_back(mindex() + 1);
		}

		_current += match.length(0);
	}

	{ // exchange block ends
		if (!Communication::allGatherUnknownSize(blockends)) {
			ESINFO(ERROR) << "ESPRESO internal error: exchange Workbench block ends.";
		}
	}

	{ // find coordinate blocks ends
		for (size_t i = 0; i < _coordinates.size(); i++) {
			// NBLOCK, NUMFIELD, Solkey
			std::string command = getLine(_coordinates[i].header);
			std::string datadesc = getLine(_coordinates[i].header + command.size() + 2);

			std::vector<std::string> nblock = Parser::split(command, ",");
			std::vector<std::string> format = Parser::split(datadesc.substr(1, datadesc.size() - 2), ",");
			std::vector<std::string> index  = Parser::split(format[0], "i");
			std::vector<std::string> coords = Parser::split(format[1], "e.");
			int NUMFIELD = std::stoi(nblock[1]);
			int indexsize = std::stoi(index[1]), datacount = std::stoi(coords[0]), datasize = std::stoi(coords[1]);
			if (nblock.size() == 3) {
				ESINFO(ERROR) << "Workbench parse error: not supported NBLOCK format (Solkey).";
			}

			_coordinates[i].indexsize = indexsize;
			_coordinates[i].datacount = datacount;
			_coordinates[i].datasize = datasize;
			_coordinates[i].sIndex = _coordinates[i].header + command.size() + datadesc.size() + 4;
			_coordinates[i].eIndex = *std::lower_bound(blockends.begin(), blockends.end(), _coordinates[i].header);
			fillranks(_coordinates[i]);
		}

		if (!Communication::allGatherUnknownSize(_coordinates)) {
			ESINFO(ERROR) << "ESPRESO internal error: exchange Workbench coordinates blocks.";
		}
	}

	{ // find elements blocks ends
		for (size_t i = 0; i < _elements.size(); i++) {
			// EBLOCK, NUM_NODES, Solkey
			std::string command = getLine(_elements[i].header);
			std::string datadesc = getLine(_elements[i].header + command.size() + 2);

			std::vector<std::string> eblock = Parser::split(command, ",");
			std::vector<std::string> format = Parser::split(datadesc.substr(1, datadesc.size() - 2), "i");
			int NUM_NODES = std::stoi(eblock[1]);
			int count = std::stoi(format[0]), size = std::stoi(format[1]);

			_elements[i].datacount = count;
			_elements[i].datasize = size;
			_elements[i].sIndex = _elements[i].header + command.size() + datadesc.size() + 4;
			_elements[i].eIndex = *std::lower_bound(blockends.begin(), blockends.end(), _elements[i].header);
			fillranks(_elements[i]);
		}

		if (!Communication::allGatherUnknownSize(_elements)) {
			ESINFO(ERROR) << "ESPRESO internal error: exchange Workbench coordinates blocks.";
		}
	}

	{ // parse CMBLOCKs
		for (size_t i = 0; i < _cmblocks.size(); i++) {
			// CMBLOCK, Cname, Entity, NUMITEMS
			std::string command = getLine(_cmblocks[i].header);
			std::string datadesc = getLine(_cmblocks[i].header + command.size() + 2);

			std::vector<std::string> cmblock = Parser::split(command, ",");
			std::vector<std::string> format = Parser::split(datadesc.substr(1, datadesc.size() - 2), "i");
			int NUMITEMS = std::stoi(cmblock[3]);
			int count = std::stoi(format[0]), size = std::stoi(format[1]);
			int lines = NUMITEMS / count;

			_cmblocks[i].datacount = count;
			_cmblocks[i].datasize = size;
			_cmblocks[i].sIndex = _cmblocks[i].header + command.size() + datadesc.size() + 4;
			_cmblocks[i].eIndex = _cmblocks[i].sIndex + NUMITEMS * size + 2 * lines + 2;
			fillranks(_cmblocks[i]);
		}

		if (!Communication::allGatherUnknownSize(_cmblocks)) {
			ESINFO(ERROR) << "ESPRESO internal error: exchange Workbench coordinates blocks.";
		}
	}

	{ // parse etypes
		for (size_t et = 0; et < etypes.size(); et++) {
			std::vector<std::string> etype = Parser::split(getLine(etypes[et]), ",");
			int index = std::stoi(etype[1]), value = std::stoi(etype[2]);
			if (_etypes.size() < index) {
				_etypes.resize(index);
			}
			_etypes[index - 1] = value;
		}

		if (!Communication::allGatherUnknownSize(_etypes)) {
			ESINFO(ERROR) << "ESPRESO internal error: exchange Workbench etype.";
		}
	}

	// TODO: make it more general (now only 20 node elements works or elements with 1 line)
	{ // fix distribution if blocks are across more processes
		for (size_t i = 0; i < _elements.size(); i++) {
			if (_elements[i].sRank != _elements[i].eRank && _elements[i].datacount > 18) {
				int size1 = _elements[i].datacount * _elements[i].datasize + 2; // first line
				int size2 = (20 - 8) * _elements[i].datasize + 2; // second line (20 node element)

				for (int rank = _elements[i].sRank + 1; rank <= _elements[i].eRank; rank++) {
					if ((_dataOffset[rank] - _elements[i].sIndex) % (size1 + size2) != 0) {
						_dataOffset[rank] += size2;
						if (rank - 1 == environment->MPIrank) {
							_end += size2;
						}
						if (rank == environment->MPIrank) {
							_begin += size2;
							_current = _begin;
						}
					}
				}
			}
		}
	}
}


