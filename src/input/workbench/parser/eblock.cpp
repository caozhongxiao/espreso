
#include "eblock.h"
#include "et.h"
#include "../../../basis/containers/tarray.h"
#include "../../../basis/utilities/parser.h"
#include "../../../basis/utilities/communication.h"
#include "../../../basis/utilities/utils.h"

#include "../../../mesh/elements/element.h"
#include "../../input.h"

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
	case 5:
		NDSEL = command[4].size() ? std::stol(command[4]) : -1;
	case 4:
		NDMAX = command[3].size() ? std::stol(command[3]) : -1;
	case 3:
		NUM_NODES = std::stoi(command[1]);
		Solkey = command[2].size();
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

void EBlock::fixOffsets(std::vector<size_t> &dataOffsets)
{
	if (Solkey) {
		eslocal nodes;
		if (fRank == environment->MPIrank) {
			nodes = std::stoi(Parser::split(Parser::getLine(begin + first - offset), " ")[9]);
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
		if (nodes < 8) {
			valueSize = valueSize - 8 + nodes;
			lineSize = valueSize * valueLength + lineEndSize;
			MPI_Bcast(&valueSize, sizeof(eslocal), MPI_BYTE, fRank, environment->MPICommunicator);
			MPI_Bcast(&lineSize, sizeof(eslocal), MPI_BYTE, fRank, environment->MPICommunicator);
			size1 = elementSize = lineSize;
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
		MPI_Bcast(&valueSize, sizeof(eslocal), MPI_BYTE, fRank, environment->MPICommunicator);
		MPI_Bcast(&lineSize, sizeof(eslocal), MPI_BYTE, fRank, environment->MPICommunicator);
		elementSize = lineSize;
	}
}

bool EBlock::readData(const std::vector<ET> &et, PlainMeshData &mesh)
{
	if (Solkey) {
		return solid(et, mesh);
	} else {
		return boundary(et, mesh);
	}
}

bool EBlock::solid(const std::vector<ET> &et, PlainMeshData &mesh)
{
	size_t threads = environment->OMP_NUM_THREADS;

	const char *first = getFirst(), *last = getLast();
	size_t size = (last - first) / elementSize;

	std::vector<size_t> tdistribution = tarray<size_t>::distribute(threads, size);

	std::vector<std::vector<eslocal> > tesize(threads), tnodes(threads), tIDs(threads);
	std::vector<std::vector<int> > ttype(threads), tbody(threads), tmat(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<char> value(valueLength + 1);
		std::vector<eslocal> nindices(20);
		int nnodes;

		for (auto element = first + elementSize * tdistribution[t]; element < first + elementSize * tdistribution[t + 1];) {
			tbody[t].push_back(0);
			tmat[t].push_back(atoi(element) - 1); element += valueLength; // material
			ttype[t].push_back(atoi(element) - 1); element += valueLength; // etype
			element += valueLength; // real constant
			element += valueLength; // section ID
			element += valueLength; // element coordinate system
			element += valueLength; // birth / death
			element += valueLength; // solid model reference number
			element += valueLength; // element shape flag
			nnodes = atoi(element); element += valueLength; // number of nodes
			element += valueLength; // not used
			tIDs[t].push_back(atoi(element) - 1); element += valueLength; // element ID

			for (int i = 0; i < 8 && i < nnodes; i++) {
				memcpy(value.data(), element, valueLength);
				element += valueLength; // element ID
				nindices[i] = atol(value.data()) - 1;
			}
			element += lineEndSize;


			auto readNextNodes = [&] () {
				for (int i = 0; i < nnodes - 8; i++) {
					memcpy(value.data(), element, valueLength);
					element += valueLength;
					nindices[i + 8] = atol(value.data()) - 1;
				}
				element += lineEndSize;
			};

			switch (et[ttype[t].back()].etype()) {
			case ET::ETYPE::D2SOLID_4NODES:
				if (nindices[2] == nindices[3]) { // triangle3
					tesize[t].push_back(3);
					ttype[t].back() = (eslocal)Element::CODE::TRIANGLE3;
					tnodes[t].insert(tnodes[t].end(), nindices.begin(), nindices.begin() + 3);
				} else { // square4
					tesize[t].push_back(4);
					ttype[t].back() = (eslocal)Element::CODE::SQUARE4;
					tnodes[t].insert(tnodes[t].end(), nindices.begin(), nindices.begin() + 4);
				}
				break;
			case ET::ETYPE::D2SOLID_6NODES:
				tesize[t].push_back(6);
				ttype[t].back() = (eslocal)Element::CODE::TRIANGLE6;
				tnodes[t].insert(tnodes[t].end(), nindices.begin(), nindices.begin() + 6);
				break;
			case ET::ETYPE::D2SOLID_8NODES:
				if (nindices[2] == nindices[3]) { // triangle6
					tesize[t].push_back(6);
					ttype[t].back() = (eslocal)Element::CODE::TRIANGLE6;
					tnodes[t].insert(tnodes[t].end(), nindices.begin(), nindices.begin() + 3);
					tnodes[t].insert(tnodes[t].end(), nindices.begin() + 4, nindices.begin() + 6);
					tnodes[t].push_back(nindices[7]);
				} else { // square8
					tesize[t].push_back(8);
					ttype[t].back() = (eslocal)Element::CODE::SQUARE8;
					tnodes[t].insert(tnodes[t].end(), nindices.begin(), nindices.begin() + 8);
				}
				break;

			case ET::ETYPE::D3SOLID_4NODES:
				tesize[t].push_back(4);
				ttype[t].back() = (eslocal)Element::CODE::TETRA4;
				tnodes[t].insert(tnodes[t].end(), nindices.begin(), nindices.begin() + 4);
				break;
			case ET::ETYPE::D3SOLID_8NODES:
				if (nindices[2] == nindices[3]) {
					if (nindices[4] == nindices[5]) { // tetra4
						tesize[t].push_back(4);
						ttype[t].back() = (eslocal)Element::CODE::TETRA4;
						tnodes[t].insert(tnodes[t].end(), nindices.begin(), nindices.begin() + 3);
						tnodes[t].push_back(nindices[4]);
					} else { // prisma6
						tesize[t].push_back(6);
						ttype[t].back() = (eslocal)Element::CODE::PRISMA6;
						tnodes[t].insert(tnodes[t].end(), nindices.begin(), nindices.begin() + 3);
						tnodes[t].insert(tnodes[t].end(), nindices.begin() + 4, nindices.begin() + 7);
					}
				} else {
					if (nindices[4] == nindices[5]) { // pyramid5
						tesize[t].push_back(5);
						ttype[t].back() = (eslocal)Element::CODE::PYRAMID5;
						tnodes[t].insert(tnodes[t].end(), nindices.begin(), nindices.begin() + 5);
					} else { // hexa8
						tesize[t].push_back(8);
						ttype[t].back() = (eslocal)Element::CODE::HEXA8;
						tnodes[t].insert(tnodes[t].end(), nindices.begin(), nindices.begin() + 8);
					}
				}
				break;
			case ET::ETYPE::D3SOLID_10NODES:
				readNextNodes();
				tesize[t].push_back(10);
				ttype[t].back() = (eslocal)Element::CODE::TETRA10;
				tnodes[t].insert(tnodes[t].end(), nindices.begin(), nindices.begin() + 10);
				break;
			case ET::ETYPE::D3SOLID_20NODES:
				readNextNodes();
				if (nindices[2] == nindices[3]) {
					if (nindices[4] == nindices[5]) { // tetra10
						tesize[t].push_back(10);
						ttype[t].back() = (eslocal)Element::CODE::TETRA10;
						tnodes[t].insert(tnodes[t].end(), nindices.begin(), nindices.begin() + 3);
						tnodes[t].push_back(nindices[4]);

						tnodes[t].insert(tnodes[t].end(), nindices.begin() + 8, nindices.begin() + 10);
						tnodes[t].push_back(nindices[11]);
						tnodes[t].insert(tnodes[t].end(), nindices.begin() + 16, nindices.begin() + 19);
					} else { // prisma15
						tesize[t].push_back(15);
						ttype[t].back() = (eslocal)Element::CODE::PRISMA15;
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
						ttype[t].back() = (eslocal)Element::CODE::PYRAMID13;
						tnodes[t].insert(tnodes[t].end(), nindices.begin(), nindices.begin() + 5);
						tnodes[t].insert(tnodes[t].end(), nindices.begin() + 8, nindices.begin() + 12);
						tnodes[t].insert(tnodes[t].end(), nindices.begin() + 16, nindices.begin() + 20);
					} else { // hexa20
						tesize[t].push_back(20);
						ttype[t].back() = (eslocal)Element::CODE::HEXA20;
						tnodes[t].insert(tnodes[t].end(), nindices.begin(), nindices.begin() + 20);
					}
				}
				break;
			default:
				ESINFO(ERROR) << "ESPRESO Workbench parser: not implemented parsing of etype: " << ttype[t].back() << " = type " << et[ttype[t].back()].type;
			}
		}
	}

	for (size_t t = 0; t < threads; t++) {
		mesh.eIDs.insert(mesh.eIDs.end(), tIDs[t].begin(), tIDs[t].end());
		mesh.esize.insert(mesh.esize.end(), tesize[t].begin(), tesize[t].end());
		mesh.enodes.insert(mesh.enodes.end(), tnodes[t].begin(), tnodes[t].end());
		mesh.etype.insert(mesh.etype.end(), ttype[t].begin(), ttype[t].end());
		mesh.body.insert(mesh.body.end(), tbody[t].begin(), tbody[t].end());
		mesh.material.insert(mesh.material.end(), tmat[t].begin(), tmat[t].end());
	}

	return true;
}

bool EBlock::boundary(const std::vector<ET> &et, PlainMeshData &mesh)
{
	size_t threads = environment->OMP_NUM_THREADS;

	const char *first = getFirst(), *last = getLast();
	eslocal size = (last - first) / elementSize;

	std::vector<size_t> tdistribution = tarray<eslocal>::distribute(threads, size);

	std::vector<std::vector<eslocal> > tesize(threads), tnodes(threads), tIDs(threads);
	std::vector<std::vector<int> > ttype(threads), tbody(threads), tmat(threads);
	int nodes = valueSize - 5;

	if (nodes != 4 && nodes != 8 && nodes != 2 && nodes != 3) {
		ESINFO(ERROR) << "ESPRESO Workbench parser: uknown format of EBLOCK. Nodes = " << nodes;
	}

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<char> value(valueLength + 1);
		std::vector<eslocal> nindices(20);

		for (auto element = first + elementSize * tdistribution[t]; element < first + elementSize * tdistribution[t + 1];) {
			tIDs[t].push_back(atoi(element) - 1); element += valueLength; // element ID
			tbody[t].push_back(0);
			ttype[t].push_back(0);
			element += valueLength; // section ID
			element += valueLength; // real constant
			tmat[t].push_back(atoi(element) - 1); element += valueLength; // material
			element += valueLength; // element coordinate system

			for (int i = 0; i < nodes; i++) {
				memcpy(value.data(), element, valueLength);
				element += valueLength; // element node
				nindices[i] = atol(value.data()) - 1;
			}
			element += lineEndSize;

			if (nodes == 2) { // line2
				tesize[t].push_back(2);
				ttype[t].back() = (eslocal)Element::CODE::LINE2;
				tnodes[t].insert(tnodes[t].end(), nindices.begin(), nindices.begin() + 2);
			}

			if (nodes == 3) { // line3
				tesize[t].push_back(3);
				ttype[t].back() = (eslocal)Element::CODE::LINE3;
				tnodes[t].insert(tnodes[t].end(), nindices.begin(), nindices.begin() + 3);
			}

			if (nodes == 4) {
				if (nindices[2] == nindices[3]) { // triangle3
					tesize[t].push_back(3);
					ttype[t].back() = (eslocal)Element::CODE::TRIANGLE3;
					tnodes[t].insert(tnodes[t].end(), nindices.begin(), nindices.begin() + 3);
				} else { // square4
					tesize[t].push_back(4);
					ttype[t].back() = (eslocal)Element::CODE::SQUARE4;
					tnodes[t].insert(tnodes[t].end(), nindices.begin(), nindices.begin() + 4);
				}
			}
			if (nodes == 8) {
				if (nindices[2] == nindices[3]) { // triangle6
					tesize[t].push_back(6);
					ttype[t].back() = (eslocal)Element::CODE::TRIANGLE6;
					tnodes[t].insert(tnodes[t].end(), nindices.begin(), nindices.begin() + 3);
					tnodes[t].insert(tnodes[t].end(), nindices.begin() + 4, nindices.begin() + 6);
					tnodes[t].push_back(nindices[7]);
				} else { // square8
					tesize[t].push_back(8);
					ttype[t].back() = (eslocal)Element::CODE::SQUARE8;
					tnodes[t].insert(tnodes[t].end(), nindices.begin(), nindices.begin() + 8);
				}
			}
		}
	}

	for (size_t t = 0; t < threads; t++) {
		mesh.eIDs.insert(mesh.eIDs.end(), tIDs[t].begin(), tIDs[t].end());
		mesh.esize.insert(mesh.esize.end(), tesize[t].begin(), tesize[t].end());
		mesh.enodes.insert(mesh.enodes.end(), tnodes[t].begin(), tnodes[t].end());
		mesh.etype.insert(mesh.etype.end(), ttype[t].begin(), ttype[t].end());
		mesh.body.insert(mesh.body.end(), tbody[t].begin(), tbody[t].end());
		mesh.material.insert(mesh.material.end(), tmat[t].begin(), tmat[t].end());
	}

	return true;
}


