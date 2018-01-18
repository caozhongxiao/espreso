
#include "workbench.h"

#include "../../basis/containers/point.h"
#include "../../basis/containers/serializededata.h"
#include "../../basis/logging/logging.h"
#include "../../basis/utilities/parser.h"
#include "../../basis/utilities/communication.h"
#include "../../basis/utilities/utils.h"

#include "../../config/ecf/ecf.h"

#include "../../mesh/mesh.h"
#include "../../mesh/elements/element.h"
#include "../../mesh/store/nodestore.h"
#include "../../mesh/store/elementstore.h"
#include "../../mesh/store/elementsregionstore.h"
#include "../../mesh/store/boundaryregionstore.h"

#include <iostream>
#include <string>
#include <regex>
#include <cstdlib>
#include <numeric>

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

	std::vector<Point> coordinates;
	std::vector<eslocal> edist, nodes, data;

	readCoordinates(coordinates);
	readElements(edist, nodes, data);

	fillMesh(coordinates, edist, nodes, data);

	addNodeRegions();
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

			std::vector<std::string> nblock = Parser::split(command, ",", false);
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
			_coordinates[i].size = (_coordinates[i].eIndex - _coordinates[i].sIndex) / (indexsize + datacount * datasize + 2);
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

			std::vector<std::string> eblock = Parser::split(command, ",", false);
			std::vector<std::string> format = Parser::split(datadesc.substr(1, datadesc.size() - 2), "i");
			int NUM_NODES = std::stoi(eblock[1]), n;
			int count = std::stoi(format[0]), size = std::stoi(format[1]);

			if (eblock.size() == 5) {
				_elements[i].size = std::stoi(eblock[4]);
			}

			_elements[i].datacount = count;
			_elements[i].datasize = size;
			_elements[i].sIndex = _elements[i].header + command.size() + datadesc.size() + 4;
			_elements[i].eIndex = *std::lower_bound(blockends.begin(), blockends.end(), _elements[i].header);
			fillranks(_elements[i]);
			if (!StringCompare::caseInsensitiveEq(eblock[2], "solid")) {
				_boundary.push_back(_elements[i]);
				_elements.erase(_elements.begin() + i);
			}
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

			std::vector<std::string> cmblock = Parser::split(command, ",", false);
			std::vector<std::string> format = Parser::split(datadesc.substr(1, datadesc.size() - 2), "i");
			int NUMITEMS = std::stoi(cmblock[3]);
			int count = std::stoi(format[0]), size = std::stoi(format[1]);
			int lines = NUMITEMS / count;

			_cmblocks[i].size = NUMITEMS;
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
			std::vector<std::string> etype = Parser::split(getLine(etypes[et]), ",", false);
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
			if (_elements[i].sRank == environment->MPIrank) {
				_elements[i].nodes = std::stoi(Parser::split(getLine(_elements[i].sIndex), " ")[9]);
			}
			MPI_Bcast(&_elements[i].nodes, sizeof(eslocal), MPI_BYTE, _elements[i].sRank, environment->MPICommunicator);

			int size1 = _elements[i].datacount * _elements[i].datasize + 2; // first line
			int size2 = 0;
			if (_elements[i].nodes > 8) {
				size2 = (_elements[i].nodes - 8) * _elements[i].datasize + 2; // second line (20 node element)
			}
			if (_elements[i].sRank != _elements[i].eRank && _elements[i].nodes > 8) {
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
			if (_elements[i].size == -1) {
				_elements[i].size = (_elements[i].eIndex - _elements[i].sIndex) / (size1 + size2);
			}
		}
	}
}

void WorkbenchLoader::interval(DataInterval &interval, char* &first, char* &last)
{
	first = _begin;
	last = _end;
	if (interval.sRank == environment->MPIrank) {
		first = _begin + interval.sIndex - _dataOffset[environment->MPIrank];
	}
	if (interval.eRank == environment->MPIrank) {
		last = _begin + interval.eIndex - _dataOffset[environment->MPIrank];
	}
}

void WorkbenchLoader::readCoordinates(std::vector<Point> &coordinates)
{
	if (_coordinates.size() != 1) {
		ESINFO(GLOBAL_ERROR) << "Workbench parse error: implement parsing of more nblocks.";
	}

	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<size_t> ndistribution = tarray<eslocal>::distribute(environment->MPIsize, _coordinates[0].size);

	std::vector<int> targets;
	std::vector<std::vector<Point> > sBuffer, rBuffer;

	std::vector<std::vector<Point> > nodes(threads);
	eslocal offset, size;
	for (size_t i = 0; i < _coordinates.size(); i++) {
		if (_coordinates[i].sRank <= environment->MPIrank && environment->MPIrank <= _coordinates[i].eRank) {
			char *first, *last;
			interval(_coordinates[i], first, last);

			int linesize = _coordinates[i].indexsize + _coordinates[i].datacount * _coordinates[i].datasize + 2;
			offset = (first - _begin + _dataOffset[environment->MPIrank] - _coordinates[i].sIndex) / linesize;
			size = (last - first) / linesize;

			std::vector<size_t> tdistribution = tarray<eslocal>::distribute(threads, size);

			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				nodes[t].reserve(tdistribution[t + 1] - tdistribution[t]);
				double x, y, z;
				if (_coordinates[i].datacount == 3) {
					for (auto data = first + linesize * tdistribution[t]; data != first + linesize * tdistribution[t + 1];) {
						data += _coordinates[i].indexsize; // skip index

						x = atof(data); data += _coordinates[i].datasize;
						y = atof(data); data += _coordinates[i].datasize;
						z = atof(data); data += _coordinates[i].datasize;

						data += 2;

						nodes[t].push_back(Point(x, y, z));
					}
				}
				if (_coordinates[i].datacount == 6) {
					for (auto data = first + linesize * tdistribution[t]; data != first + linesize * tdistribution[t + 1];) {
						data += _coordinates[i].indexsize; // skip index

						x = atof(data); data += _coordinates[i].datasize;
						y = atof(data); data += _coordinates[i].datasize;
						z = atof(data); data += _coordinates[i].datasize;

						data += _coordinates[i].datasize;
						data += _coordinates[i].datasize;
						data += _coordinates[i].datasize;

						data += 2;

						nodes[t].push_back(Point(x, y, z));
					}
				}
			}
		}
	}

	eslocal nextindex = offset, t = 0, tindex = 0;
	for (int rank = 0; rank < environment->MPIsize && nextindex < offset + size; rank++) {
		if (ndistribution[rank] <= nextindex) {
			sBuffer.push_back({});
			targets.push_back(rank);
			while (nextindex < ndistribution[rank + 1]) {
				if (nodes[t].size() - tindex <= ndistribution[rank + 1] - nextindex) {
					sBuffer.back().insert(sBuffer.back().end(), nodes[t].begin() + tindex, nodes[t].end());
					nextindex += nodes[t++].size() - tindex;
					tindex = 0;
				} else {
					eslocal size = nodes[t].size() - tindex;
					sBuffer.back().insert(sBuffer.back().end(), nodes[t].begin() + tindex, nodes[t].begin() + tindex + size);
					nextindex += size;
					tindex += size;
				}
			}
		}
	}

	if (!Communication::sendVariousTargets(sBuffer, rBuffer, targets)) {
		ESINFO(ERROR) << "ESPRESO internal error: distribute Workbench example nodes.";
	}

	for (size_t r = 0; r < rBuffer.size(); r++) {
		coordinates.insert(coordinates.end(), rBuffer[r].begin(), rBuffer[r].end());
	}
}

void WorkbenchLoader::readElements(std::vector<eslocal> &edist, std::vector<eslocal> &nodes, std::vector<eslocal> &data)
{
	size_t threads = environment->OMP_NUM_THREADS;

	eslocal esize = 0;
	for (size_t eb = 0; eb < _elements.size(); eb++) {
		esize += _elements[eb].size;
	}

	std::vector<int> targets;
	std::vector<std::vector<eslocal> > sBuffer, rBuffer;

	std::vector<size_t> edistribution = tarray<eslocal>::distribute(environment->MPIsize, esize);

	std::vector<std::vector<eslocal> > elements(threads);
	for (size_t eb = 0, eboffset = 0; eb < _elements.size(); eboffset += _elements[eb++].size) {
		int linesize = _elements[eb].datacount * _elements[eb].datasize + 2;
		if (_elements[eb].nodes > 8) {
			linesize += (_elements[eb].nodes - 8) * _elements[eb].datasize + 2;
		}

		if (_elements[eb].sRank <= environment->MPIrank && environment->MPIrank <= _elements[eb].eRank) {
			char *first, *last;
			interval(_elements[eb], first, last);

			eslocal offset = eboffset + (first - _begin + _dataOffset[environment->MPIrank] - _elements[eb].sIndex) / linesize;
			eslocal size = (last - first) / linesize;

			int ftarget = std::lower_bound(edistribution.begin(), edistribution.end(), offset + 1) - edistribution.begin() - 1;
			int ltarget = std::lower_bound(edistribution.begin(), edistribution.end(), offset + size + 1) - edistribution.begin() - 1;

			for (int target = ftarget; target < ltarget; target++) {
				eslocal tsize = edistribution[target + 1] - offset;
				if (tsize > size) {
					tsize = size;
				}

				std::vector<size_t> tdistribution = tarray<eslocal>::distribute(threads, tsize);

				#pragma omp parallel for
				for (size_t t = 0; t < threads; t++) {
					std::vector<eslocal> nodes(20);
					std::vector<char> nindex(_elements[eb].datasize + 1);
					int nnodes;
					int datasize = 3 + _elements[eb].nodes; // region, mat, code, nodes
					elements[t].reserve(elements[t].size() + (tdistribution[t + 1] - tdistribution[t]) * datasize);

					for (auto edata = first + linesize * tdistribution[t]; edata != first + linesize * tdistribution[t + 1];) {
						elements[t].push_back(eb);
						elements[t].push_back(atoi(edata) - 1); edata += _elements[eb].datasize; // material
						edata += _elements[eb].datasize; // etype
						edata += _elements[eb].datasize; // real constant
						edata += _elements[eb].datasize; // section ID
						edata += _elements[eb].datasize; // element coordinate system
						edata += _elements[eb].datasize; // birth / death
						edata += _elements[eb].datasize; // solid model reference number
						edata += _elements[eb].datasize; // element shape flag
						nnodes = atoi(edata); edata += _elements[eb].datasize; // number of nodes
						edata += _elements[eb].datasize; // not used
						edata += _elements[eb].datasize; // element ID

						for (int i = 0; i < 8; i++) {
							memcpy(nindex.data(), edata, _elements[eb].datasize);
							nodes[i] = atol(nindex.data()) - 1; edata += _elements[eb].datasize; // element ID
						}
						edata += 2;

						if (nnodes != _elements[eb].nodes) {
							ESINFO(ERROR) << "Workbench parse error: elements in eblock have various number of nodes.";
						}
						if (nnodes > 8) {
							for (int i = 0; i < _elements[eb].nodes - 8; i++) {
								memcpy(nindex.data(), edata, _elements[eb].datasize);
								nodes[i + 8] = atol(nindex.data()) - 1; edata += _elements[eb].datasize; // element ID
							}
							edata += 2;
						}

						if (nnodes == 20) {
							if (nodes[2] == nodes[3]) {
								if (nodes[4] == nodes[5]) { // tetra10
									elements[t].push_back((eslocal)Element::CODE::TETRA10);
									elements[t].insert(elements[t].end(), nodes.begin(), nodes.begin() + 3);
									elements[t].push_back(nodes[4]);

									elements[t].insert(elements[t].end(), nodes.begin() + 8, nodes.begin() + 10);
									elements[t].push_back(nodes[11]);
									elements[t].insert(elements[t].end(), nodes.begin() + 16, nodes.begin() + 19);
								} else { // prisma15
									elements[t].push_back((eslocal)Element::CODE::PRISMA15);
									elements[t].insert(elements[t].end(), nodes.begin(), nodes.begin() + 3);
									elements[t].insert(elements[t].end(), nodes.begin() + 4, nodes.begin() + 7);

									elements[t].insert(elements[t].end(), nodes.begin() + 8, nodes.begin() + 10);
									elements[t].push_back(nodes[11]);
									elements[t].insert(elements[t].end(), nodes.begin() + 12, nodes.begin() + 14);
									elements[t].push_back(nodes[15]);
									elements[t].insert(elements[t].end(), nodes.begin() + 16, nodes.begin() + 19);
								}
							} else {
								if (nodes[4] == nodes[5]) { // pyramid13
									elements[t].push_back((eslocal)Element::CODE::PYRAMID13);
									elements[t].insert(elements[t].end(), nodes.begin(), nodes.begin() + 5);
									elements[t].insert(elements[t].end(), nodes.begin() + 8, nodes.begin() + 12);
									elements[t].insert(elements[t].end(), nodes.begin() + 16, nodes.begin() + 20);
								} else { // hexa20
									elements[t].push_back((eslocal)Element::CODE::HEXA20);
									elements[t].insert(elements[t].end(), nodes.begin(), nodes.begin() + 20);
								}
							}
						}
						if (nnodes == 8) {
							if (nodes[2] == nodes[3]) {
								if (nodes[4] == nodes[5]) { // tetra4
									elements[t].push_back((eslocal)Element::CODE::TETRA4);
									elements[t].insert(elements[t].end(), nodes.begin(), nodes.begin() + 3);
									elements[t].push_back(nodes[4]);
								} else { // prisma6
									elements[t].push_back((eslocal)Element::CODE::PRISMA6);
									elements[t].insert(elements[t].end(), nodes.begin(), nodes.begin() + 3);
									elements[t].insert(elements[t].end(), nodes.begin() + 4, nodes.begin() + 7);
								}
							} else {
								if (nodes[4] == nodes[5]) { // pyramid5
									elements[t].push_back((eslocal)Element::CODE::PYRAMID5);
									elements[t].insert(elements[t].end(), nodes.begin(), nodes.begin() + 5);
								} else { // hexa8
									elements[t].push_back((eslocal)Element::CODE::HEXA8);
									elements[t].insert(elements[t].end(), nodes.begin(), nodes.begin() + 8);
								}
							}
						}
						if (nnodes == 10) {
							elements[t].push_back((eslocal)Element::CODE::TETRA10);
							elements[t].insert(elements[t].end(), nodes.begin(), nodes.begin() + 10);
						}
					}
				}

				first += tsize * linesize;
				offset += tsize;
				size -= tsize;
				targets.push_back(target);
				sBuffer.push_back({});
				for (size_t t = 0; t < threads; t++) {
					sBuffer.back().insert(sBuffer.back().end(), elements[t].begin(), elements[t].end());
					elements[t].clear();
				}
			}
		}

		if (!Communication::sendVariousTargets(sBuffer, rBuffer, targets)) {
			ESINFO(ERROR) << "ESPRESO internal error: distribute Workbench example nodes.";
		}

		edist.push_back(0);
		for (size_t t = 0; t < rBuffer.size(); t++) {
			size_t i = 0;
			while (i < rBuffer[t].size()) {
				data.push_back(rBuffer[t][i++]); // block
				data.push_back(rBuffer[t][i++]); // material
				data.push_back(rBuffer[t][i]); // code
				edist.push_back(edist.back() + _mesh._eclasses[0][rBuffer[t][i]].nodes);
				for (size_t n = 0; n < _mesh._eclasses[0][rBuffer[t][i]].nodes; n++) {
					nodes.push_back(rBuffer[t][i + n + 1]);
				}
				i += _mesh._eclasses[0][rBuffer[t][i]].nodes + 1;
			}
		}
	}
}

void WorkbenchLoader::addNodeRegions()
{
	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<int> targets;
	std::vector<std::vector<eslocal> > sBuffer, rBuffer;

	std::vector<std::vector<eslocal> > tnodes(threads);
	eslocal offset, size;
	for (size_t i = 0; i < _cmblocks.size(); i++) {
		if (_cmblocks[i].sRank <= environment->MPIrank && environment->MPIrank <= _cmblocks[i].eRank) {
			char *first, *last;
			interval(_cmblocks[i], first, last);

			int linesize = _cmblocks[i].datacount * _cmblocks[i].datasize + 2;
			offset = (first - _begin + _dataOffset[environment->MPIrank] - _cmblocks[i].sIndex) / linesize;
			size = (last - first) / linesize;

			std::vector<size_t> tdistribution = tarray<eslocal>::distribute(threads, size);

			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				std::vector<char> nindex(_cmblocks[i].datasize + 1);
				tnodes[t].reserve((tdistribution[t + 1] - tdistribution[t]) * _cmblocks[i].datacount);
				for (auto edata = first + linesize * tdistribution[t]; edata != first + linesize * tdistribution[t + 1];) {
					for (eslocal n = 0; n < _cmblocks[i].datacount; ++n) {
						memcpy(nindex.data(), edata, _cmblocks[i].datasize);
						tnodes[t].push_back(atol(nindex.data()) - 1); edata += _cmblocks[i].datasize;
					}
					edata += 2;
				}
			}
			first += linesize * tdistribution[threads];

			while (first != last - 2) {
				std::vector<char> nindex(_cmblocks[i].datasize + 1);
				memcpy(nindex.data(), first, _cmblocks[i].datasize);
				tnodes.back().push_back(atol(nindex.data()) - 1); first += _cmblocks[i].datasize;
			}
		}
	}

	serializededata<eslocal, eslocal>::balance(1, tnodes);
	_mesh.boundaryRegions.push_back(new BoundaryRegionStore("FIX", _mesh._eclasses));
	_mesh.boundaryRegions.back()->nodes = new serializededata<eslocal, eslocal>(1, tnodes);
}

void WorkbenchLoader::fillMesh(std::vector<Point> &coordinates, std::vector<eslocal> &edist, std::vector<eslocal> &nodes, std::vector<eslocal> &data)
{
	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<std::vector<eslocal> > nIDs(threads), nRanks(threads);
	std::vector<std::vector<Point> > tcoordinates(threads);

	std::vector<std::vector<eslocal> > tedist(threads);
	std::vector<std::vector<eslocal> > tnodes(threads);
	std::vector<std::vector<eslocal> > eIDs(threads), eMat(threads), eBody(threads), rData(threads);
	std::vector<std::vector<Element*> > epointers(threads);

	std::vector<size_t> cdistribution = tarray<Point>::distribute(threads, coordinates.size());
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		tcoordinates[t].insert(tcoordinates[t].end(), coordinates.begin() + cdistribution[t], coordinates.begin() + cdistribution[t + 1]);
		nIDs[t].resize(cdistribution[t + 1] - cdistribution[t]);
		std::iota(nIDs[t].begin(), nIDs[t].end(), cdistribution[t]);
		nRanks[t].resize(cdistribution[t + 1] - cdistribution[t]);
	}

	std::vector<size_t> edistribution = tarray<Point>::distribute(threads, edist.size() - 1);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		if(t == 0) {
			tedist[t].insert(tedist[t].end(), edist.begin() + edistribution[t], edist.begin() + edistribution[t + 1] + 1);
		} else {
			tedist[t].insert(tedist[t].end(), edist.begin() + edistribution[t] + 1, edist.begin() + edistribution[t + 1] + 1);
		}
		tnodes[t].insert(tnodes[t].end(), nodes.begin() + edist[edistribution[t]], nodes.begin() + edist[edistribution[t + 1]]);
		eIDs[t].resize(edistribution[t + 1] - edistribution[t]);
		std::iota(eIDs[t].begin(), eIDs[t].end(), edistribution[t]);
		epointers[t].reserve(edistribution[t + 1] - edistribution[t]);
		eMat[t].reserve(edistribution[t + 1] - edistribution[t]);
		eBody[t].resize(edistribution[t + 1] - edistribution[t]);
		for (size_t e = edistribution[t]; e < edistribution[t + 1]; ++e) {
			eMat[t].push_back(data[3 * e + 1]);
			epointers[t].push_back(&_mesh._eclasses[t][data[3 * e + 2]]);
		}
	}

	_mesh.nodes->size = coordinates.size();
	_mesh.nodes->distribution = cdistribution;
	_mesh.nodes->IDs = new serializededata<eslocal, eslocal>(1, nIDs);
	_mesh.nodes->ranks = new serializededata<eslocal, int>(1, nRanks);
	_mesh.nodes->coordinates = new serializededata<eslocal, Point>(1, tcoordinates);

	_mesh.elements->size = edist.size() - 1;
	_mesh.elements->distribution = edistribution;
	_mesh.elements->IDs = new serializededata<eslocal, eslocal>(1, eIDs);
	_mesh.elements->nodes = new serializededata<eslocal, eslocal>(tedist, tnodes);
	_mesh.elements->epointers = new serializededata<eslocal, Element*>(1, epointers);
	_mesh.elements->material = new serializededata<eslocal, int>(1, eMat);
	_mesh.elements->body = new serializededata<eslocal, int>(1, eBody);


	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		rData[t].resize(edistribution[t + 1] - edistribution[t]);
		std::iota(rData[t].begin(), rData[t].end(), edistribution[t]);
	}
	_mesh.elementsRegions.push_back(new ElementsRegionStore("ALL_ELEMENTS"));
	_mesh.elementsRegions.back()->elements = new serializededata<eslocal, eslocal>(1, rData);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		rData[t].resize(cdistribution[t + 1] - cdistribution[t]);
		std::iota(rData[t].begin(), rData[t].end(), cdistribution[t]);
	}
	_mesh.boundaryRegions.push_back(new BoundaryRegionStore("ALL_NODES", _mesh._eclasses));
	_mesh.boundaryRegions.back()->nodes = new serializededata<eslocal, eslocal>(1, rData);

	_mesh.neighboursWithMe.push_back(environment->MPIrank);
	std::sort(_mesh.neighboursWithMe.begin(), _mesh.neighboursWithMe.end());
}


