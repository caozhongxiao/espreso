
#include "eblock.h"
#include "../../loader.h"

#include "../../../basis/containers/tarray.h"
#include "../../../basis/utilities/parser.h"
#include "../../../basis/utilities/communication.h"
#include "../../../basis/utilities/utils.h"

#include "../../../mesh/elements/element.h"

using namespace espreso;

size_t EBlock::size = 6;
const char* EBlock::upper = "EBLOCK";
const char* EBlock::lower = "eblock";

EBlock::EBlock()
: NUM_NODES(-1), Solkey(0), NDMAX(-1), NDSEL(-1),
  lineSize(-1), elementSize(-1), lineEndSize(-1),
  valueSize(-1), valueLength(-1)
{

}

EBlock& EBlock::parse(const char* begin)
{
	auto error = [&] (std::string &line) {
		ESINFO(ERROR) << "Workbench parse error: unknown format of EBLOCK: " << line;
	};

	std::string commandLine = Parser::getLine(begin);
	std::string descriptionLine = Parser::getLine(begin + commandLine.size());

	lineEndSize = (*(commandLine.end() - 2) == '\r') ? 2 : 1;

	std::vector<std::string> command = Parser::split(commandLine, ",", false);
	std::vector<std::string> description = Parser::split(descriptionLine.substr(1, descriptionLine.size() - 2 - lineEndSize), "i");
	if (description.size() != 2) {
		error(descriptionLine);
	}

	switch (command.size()) {
	case 3:
		NUM_NODES = std::stoi(command[1]);
		Solkey = command[2].size();
		break;
	case 5:
		NUM_NODES = std::stoi(command[1]);
		Solkey = command[2].size();
		NDMAX = command[3].size() ? std::stol(command[3]) : -1;
		NDSEL = command[4].size() ? std::stol(command[4]) : -1;
		break;
	default:
		error(commandLine);
	}

	valueSize = std::stoi(description[0]);
	valueLength = std::stoi(description[1]);

	elementSize = lineSize = valueSize * valueLength + lineEndSize;

	WorkbenchParser::fillIndices(begin, begin + commandLine.size() + descriptionLine.size());
	return *this;
}

void EBlock::fixOffsets(std::vector<eslocal> &dataOffsets)
{
	if (Solkey) {
		eslocal nodes;
		if (fRank == environment->MPIrank) {
			nodes = std::stoi(Parser::split(Parser::getLine(begin + first - dataOffsets[environment->MPIrank]), " ")[9]);
		}
		MPI_Bcast(&nodes, sizeof(eslocal), MPI_BYTE, fRank, environment->MPICommunicator);

		int size1 = lineSize; // first line
		int size2 = 0;
		if (nodes > 8) {
			size2 = (nodes - 8) * valueLength + lineEndSize;
		}
		if (fRank != lRank && nodes > 8) {
			for (int rank = fRank + 1; rank <= lRank; rank++) {
				if ((dataOffsets[rank] - first) % (size1 + size2) != 0) {
					dataOffsets[rank] += size2;
					if (rank - 1 == environment->MPIrank) {
						end += size2;
					}
					if (rank == environment->MPIrank) {
						begin += size2;
						offset += size2;
					}
				}
			}
		}
		if (NDSEL == -1) {
			NDSEL = (first - last) / (size1 + size2);
		}
		elementSize = size1 + size2;
	} else {
		if (fRank == environment->MPIrank) {
			const char *first = getFirst();
			const char *next = first;
			while (*next++ != '\n');
			if ((next - first - lineEndSize) / valueLength != valueSize) {
				valueSize = (next - first - lineEndSize) / valueLength;
				lineSize = valueSize * valueLength + lineEndSize;
			}
		}
		MPI_Bcast(&lineSize, sizeof(eslocal), MPI_BYTE, fRank, environment->MPICommunicator);
		elementSize = lineSize;
	}
}

bool EBlock::readSolid(std::vector<eslocal> &edist, std::vector<eslocal> &enodes, std::vector<EData> &edata)
{
	if (!Solkey) {
		ESINFO(ERROR) << "Workbench parser internal error: EBLOCK is not solid.";
	}
	return solid(edist, enodes, edata);
}

bool EBlock::readBoundary(std::vector<eslocal> &edist, std::vector<eslocal> &enodes, std::vector<eslocal> &edata)
{
	if (Solkey) {
		ESINFO(ERROR) << "Workbench parser internal error: EBLOCK is not boundary.";
	}
	return boundary(edist, enodes, edata);
}

bool EBlock::solid(std::vector<eslocal> &esize, std::vector<eslocal> &enodes, std::vector<EData> &edata)
{
	size_t threads = environment->OMP_NUM_THREADS;

	const char *first = getFirst(), *last = getLast();
	eslocal size = (last - first) / elementSize;

	size_t offset = esize.size();

	std::vector<size_t> tdistribution = tarray<eslocal>::distribute(threads, size);

	std::vector<std::vector<eslocal> > tesize(threads), tnodes(threads);
	std::vector<std::vector<EData> > tdata(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<char> value(valueLength + 1);
		std::vector<eslocal> nindices(20);
		int nnodes;

		for (auto element = first + elementSize * tdistribution[t]; element < first + elementSize * tdistribution[t + 1];) {
			tdata[t].push_back(EData());
			tdata[t].back().body = 0;
			tdata[t].back().material = atoi(element) - 1; element += valueLength; // material
			element += valueLength; // etype
			element += valueLength; // real constant
			element += valueLength; // section ID
			element += valueLength; // element coordinate system
			element += valueLength; // birth / death
			element += valueLength; // solid model reference number
			element += valueLength; // element shape flag
			nnodes = atoi(element); element += valueLength; // number of nodes
			element += valueLength; // not used
			element += valueLength; // element ID

			for (int i = 0; i < 8; i++) {
				memcpy(value.data(), element, valueLength);
				element += valueLength; // element ID
				nindices[i] = atol(value.data()) - 1;
			}
			element += lineEndSize;

			if (nnodes > 8) {
				for (int i = 0; i < nnodes - 8; i++) {
					memcpy(value.data(), element, valueLength);
					element += valueLength;
					nindices[i + 8] = atol(value.data()) - 1;
				}
				element += lineEndSize;
			}

			if (nnodes == 20) {
				if (nindices[2] == nindices[3]) {
					if (nindices[4] == nindices[5]) { // tetra10
						tesize[t].push_back(10);
						tdata[t].back().etype = (eslocal)Element::CODE::TETRA10;
						tnodes[t].insert(tnodes[t].end(), nindices.begin(), nindices.begin() + 3);
						tnodes[t].push_back(nindices[4]);

						tnodes[t].insert(tnodes[t].end(), nindices.begin() + 8, nindices.begin() + 10);
						tnodes[t].push_back(nindices[11]);
						tnodes[t].insert(tnodes[t].end(), nindices.begin() + 16, nindices.begin() + 19);
					} else { // prisma15
						tesize[t].push_back(15);
						tdata[t].back().etype = (eslocal)Element::CODE::PRISMA15;
						tnodes[t].insert(tnodes[t].end(), nindices.begin(), nindices.begin() + 3);
						tnodes[t].insert(tnodes[t].end(), nindices.begin() + 4, nindices.begin() + 7);

						tnodes[t].insert(tnodes[t].end(), nindices.begin() + 8, nindices.begin() + 10);
						tnodes[t].push_back(nindices[11]);
						tnodes[t].insert(tnodes[t].end(), nindices.begin() + 12, nindices.begin() + 14);
						tnodes[t].push_back(nindices[15]);
						tnodes[t].insert(tnodes[t].end(), nindices.begin() + 16, nindices.begin() + 19);
					}
				} else {
					if (nindices[4] == nindices[5]) { // pyramid13
						tesize[t].push_back(13);
						tdata[t].back().etype = (eslocal)Element::CODE::PYRAMID13;
						tnodes[t].insert(tnodes[t].end(), nindices.begin(), nindices.begin() + 5);
						tnodes[t].insert(tnodes[t].end(), nindices.begin() + 8, nindices.begin() + 12);
						tnodes[t].insert(tnodes[t].end(), nindices.begin() + 16, nindices.begin() + 20);
					} else { // hexa20
						tesize[t].push_back(20);
						tdata[t].back().etype = (eslocal)Element::CODE::HEXA20;
						tnodes[t].insert(tnodes[t].end(), nindices.begin(), nindices.begin() + 20);
					}
				}
			}
			if (nnodes == 8) {
				if (nindices[2] == nindices[3]) {
					if (nindices[4] == nindices[5]) { // tetra4
						tesize[t].push_back(4);
						tdata[t].back().etype = (eslocal)Element::CODE::TETRA4;
						tnodes[t].insert(tnodes[t].end(), nindices.begin(), nindices.begin() + 3);
						tnodes[t].push_back(nindices[4]);
					} else { // prisma6
						tesize[t].push_back(6);
						tdata[t].back().etype = (eslocal)Element::CODE::PRISMA6;
						tnodes[t].insert(tnodes[t].end(), nindices.begin(), nindices.begin() + 3);
						tnodes[t].insert(tnodes[t].end(), nindices.begin() + 4, nindices.begin() + 7);
					}
				} else {
					if (nindices[4] == nindices[5]) { // pyramid5
						tesize[t].push_back(5);
						tdata[t].back().etype = (eslocal)Element::CODE::PYRAMID5;
						tnodes[t].insert(tnodes[t].end(), nindices.begin(), nindices.begin() + 5);
					} else { // hexa8
						tesize[t].push_back(8);
						tdata[t].back().etype = (eslocal)Element::CODE::HEXA8;
						tnodes[t].insert(tnodes[t].end(), nindices.begin(), nindices.begin() + 8);
					}
				}
			}
			if (nnodes == 10) {
				tesize[t].push_back(10);
				tdata[t].back().etype = (eslocal)Element::CODE::TETRA10;
				tnodes[t].insert(tnodes[t].end(), nindices.begin(), nindices.begin() + 10);
			}
		}
	}

	for (size_t t = 0; t < threads; t++) {
		esize.insert(esize.end(), tesize[t].begin(), tesize[t].end());
		enodes.insert(enodes.end(), tnodes[t].begin(), tnodes[t].end());
		edata.insert(edata.end(), tdata[t].begin(), tdata[t].end());
	}

	return true;
}

bool EBlock::boundary(std::vector<eslocal> &esize, std::vector<eslocal> &enodes, std::vector<eslocal> &edata)
{
	size_t threads = environment->OMP_NUM_THREADS;

	const char *first = getFirst(), *last = getLast();
	eslocal size = (last - first) / elementSize;

	size_t offset = esize.size();

	std::vector<size_t> tdistribution = tarray<eslocal>::distribute(threads, size);

	std::vector<std::vector<eslocal> > tesize(threads), tnodes(threads), tdata(threads);
	int nodes = valueSize - 5;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<char> value(valueLength + 1);
		std::vector<eslocal> nindices(20);

		for (auto element = first + elementSize * tdistribution[t]; element < first + elementSize * tdistribution[t + 1];) {
			element += valueLength; // element ID
			element += valueLength; // section ID
			element += valueLength; // real constant
			element += valueLength; // material
			element += valueLength; // element coordinate system

			for (int i = 0; i < nodes; i++) {
				memcpy(value.data(), element, valueLength);
				element += valueLength; // element ID
				nindices[i] = atol(value.data()) - 1;
			}
			element += lineEndSize;

			if (nodes == 4) {
				if (nindices[2] == nindices[3]) { // triangle3
					tesize[t].push_back(3);
					tdata[t].push_back((eslocal)Element::CODE::TRIANGLE3);
					tnodes[t].insert(tnodes[t].end(), nindices.begin(), nindices.begin() + 3);
				} else { // square4
					tesize[t].push_back(4);
					tdata[t].push_back((eslocal)Element::CODE::SQUARE4);
					tnodes[t].insert(tnodes[t].end(), nindices.begin(), nindices.begin() + 4);
				}
			}
			if (nodes == 8) {
				if (nindices[2] == nindices[3]) { // triangle6
					tesize[t].push_back(6);
					tdata[t].push_back((eslocal)Element::CODE::TRIANGLE6);
					tnodes[t].insert(tnodes[t].end(), nindices.begin(), nindices.begin() + 3);
					tnodes[t].insert(tnodes[t].end(), nindices.begin() + 4, nindices.begin() + 6);
					tnodes[t].push_back(nindices[7]);
				} else { // square8
					tesize[t].push_back(8);
					tdata[t].push_back((eslocal)Element::CODE::SQUARE8);
					tnodes[t].insert(tnodes[t].end(), nindices.begin(), nindices.begin() + 8);
				}
			}
		}
	}

	for (size_t t = 0; t < threads; t++) {
		esize.insert(esize.end(), tesize[t].begin(), tesize[t].end());
		enodes.insert(enodes.end(), tnodes[t].begin(), tnodes[t].end());
		edata.insert(edata.end(), tdata[t].begin(), tdata[t].end());
	}

	return true;
}


