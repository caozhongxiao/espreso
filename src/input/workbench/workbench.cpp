
#include "workbench.h"
#include "../loader.h"

#include "../../basis/containers/tarray.h"
#include "../../basis/logging/logging.h"
#include "../../basis/utilities/communication.h"

#include "../../config/ecf/ecf.h"

#include "../../mesh/mesh.h"

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
	ESINFO(OVERVIEW) << "Load ANSYS Workbench data from '" << configuration.workbench.path << "'.";
	readData();
	ESINFO(PROGRESS2) << "Workbench:: data copied from file.";
	prepareData();
	ESINFO(PROGRESS2) << "Workbench:: data prepared for reading.";

	DistributedMesh dMesh;
	parseData(dMesh);
	ESINFO(PROGRESS2) << "Workbench:: data parsed.";

	Loader::loadDistributedMesh(dMesh, mesh);
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
		const char *exchangeStart = _current;
		if (environment->MPIrank) {
			while (*_current++ != '\n');
			exchangeStart = _current;
			for (int i = 1; i < MAX_LINE_STEP; i++) {
				while (*exchangeStart++ != '\n');
			}
		}

		std::vector<std::vector<char> > sBuffer(1, std::vector<char>(const_cast<const char*>(_data.data()), exchangeStart)), rBuffer;

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

	WorkbenchParser::offset = _dataOffset[environment->MPIrank];
	WorkbenchParser::begin = _begin;
	WorkbenchParser::end = _end;
}

void WorkbenchLoader::prepareData()
{
	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<size_t> tdistribution = tarray<char>::distribute(threads, _end - _begin);
	std::vector<std::vector<NBlock> > tNBlocks(threads);
	std::vector<std::vector<EBlock> > tEBlocks(threads);
	std::vector<std::vector<CMBlock> > tCMBlocks(threads);
	std::vector<std::vector<BlockEnd> > tBlockEnds(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		const char* tbegin = _begin + tdistribution[t];
		if (t && *(tbegin - 1) != '\n') {
			while (tbegin < _end && *tbegin++ != '\n'); // start at new line
		}

		while (tbegin < _begin + tdistribution[t + 1]) {
			if (tbegin != _begin + tdistribution[t] && memcmp(tbegin - BlockEnd::unixSize, BlockEnd::unixEnd, BlockEnd::unixSize) == 0) {
				tBlockEnds[t].push_back(BlockEnd().parse(tbegin - BlockEnd::unixSize));
			}
			if (tbegin != _begin + tdistribution[t] && memcmp(tbegin - BlockEnd::winSize, BlockEnd::winEnd, BlockEnd::winSize) == 0) {
				tBlockEnds[t].push_back(BlockEnd().parse(tbegin - BlockEnd::winSize));
			}
			if (memcmp(tbegin, NBlock::upper, NBlock::size) == 0) {
				tNBlocks[t].push_back(NBlock().parse(tbegin));
			}
			if (memcmp(tbegin, NBlock::lower, NBlock::size) == 0) {
				tNBlocks[t].push_back(NBlock().parse(tbegin));
			}
			if (memcmp(tbegin, EBlock::upper, EBlock::size) == 0) {
				tEBlocks[t].push_back(EBlock().parse(tbegin));
			}
			if (memcmp(tbegin, EBlock::lower, EBlock::size) == 0) {
				tEBlocks[t].push_back(EBlock().parse(tbegin));
			}
			if (memcmp(tbegin, CMBlock::upper, CMBlock::size) == 0) {
				tCMBlocks[t].push_back(CMBlock().parse(tbegin));
			}
			if (memcmp(tbegin, CMBlock::lower, CMBlock::size) == 0) {
				tCMBlocks[t].push_back(CMBlock().parse(tbegin));
			}
			if (memcmp(tbegin, BlockEnd::nUpper, BlockEnd::nSize) == 0) {
				tBlockEnds[t].push_back(BlockEnd().parse(tbegin));
			}
			if (memcmp(tbegin, BlockEnd::nLower, BlockEnd::nSize) == 0) {
				tBlockEnds[t].push_back(BlockEnd().parse(tbegin));
			}
			while (tbegin < _end && *tbegin++ != '\n');
		}
		if (t == threads - 1) {
			if (memcmp(tbegin - BlockEnd::unixSize, BlockEnd::unixEnd, BlockEnd::unixSize) == 0) {
				tBlockEnds[t].push_back(BlockEnd().parse(tbegin - BlockEnd::unixSize));
			}
			if (memcmp(tbegin - BlockEnd::winSize, BlockEnd::winEnd, BlockEnd::winSize) == 0) {
				tBlockEnds[t].push_back(BlockEnd().parse(tbegin - BlockEnd::winSize));
			}
		}
	}

	for (size_t t = 0; t < threads; t++) {
		_NBlocks.insert(_NBlocks.end(), tNBlocks[t].begin(), tNBlocks[t].end());
		_EBlocks.insert(_EBlocks.end(), tEBlocks[t].begin(), tEBlocks[t].end());
		_CMBlocks.insert(_CMBlocks.end(), tCMBlocks[t].begin(), tCMBlocks[t].end());
		_blockEnds.insert(_blockEnds.end(), tBlockEnds[t].begin(), tBlockEnds[t].end());
	}

	if (!Communication::allGatherUnknownSize(_blockEnds)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange Workbench block ends.";
	}

	for (size_t i = 0; i < _NBlocks.size(); i++) {
		_NBlocks[i].fillDistribution(_blockEnds, _dataOffset);
	}
	for (size_t i = 0; i < _EBlocks.size(); i++) {
		_EBlocks[i].fillDistribution(_blockEnds, _dataOffset);
	}
	for (size_t i = 0; i < _CMBlocks.size(); i++) {
		_CMBlocks[i].fillDistribution(_blockEnds, _dataOffset);
	}

	if (!Communication::allGatherUnknownSize(_NBlocks)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange Workbench NBblocks.";
	}
	if (!Communication::allGatherUnknownSize(_EBlocks)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange Workbench EBlocks.";
	}
	if (!Communication::allGatherUnknownSize(_CMBlocks)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange Workbench CMBlocks.";
	}

	// fix distribution if EBlocks are across more processes and elements data have more lines
	for (size_t i = 0; i < _EBlocks.size(); i++) {
		_EBlocks[i].fixOffsets(_dataOffset);
		if (_begin != WorkbenchParser::begin) {
			_dataOffset[environment->MPIrank] += _begin - WorkbenchParser::begin;
		}
		if (_end != WorkbenchParser::end) {
			_dataOffset[environment->MPIrank + 1] += _end - WorkbenchParser::end;
		}
		_begin = WorkbenchParser::begin;
		_end = WorkbenchParser::end;
	}

	// EBlock without SOLID key are boundary blocks
	for (size_t i = _EBlocks.size() - 1; i < _EBlocks.size(); i--) {
		if (!_EBlocks[i].Solkey) {
			_BBlocks.push_back(_EBlocks[i]);
			_EBlocks.erase(_EBlocks.begin() + i);
		}
	}
}

void WorkbenchLoader::parseData(DistributedMesh &dMesh)
{
	for (size_t i = 0; i < _NBlocks.size(); i++) {
		if (!_NBlocks[i].readData(dMesh.nIDs, dMesh.coordinates)) {
			ESINFO(ERROR) << "Workbench parser: something wrong happens while read NBLOCK.";
		}
	}

	for (size_t i = 0; i < _EBlocks.size(); i++) {
		if (!_EBlocks[i].readSolid(dMesh.esize, dMesh.enodes, dMesh.edata)) {
			ESINFO(ERROR) << "Workbench parser: something wrong happens while read EBLOCK.";
		}
	}

	for (size_t i = 0; i < _BBlocks.size(); i++) {
		dMesh.bregions.push_back(MeshBRegion());
		dMesh.bregions.back().name = "BOUNDARY" + std::to_string(i + 1);
		if (!_BBlocks[i].readBoundary(dMesh.bregions.back().esize, dMesh.bregions.back().enodes, dMesh.bregions.back().etypes)) {
			ESINFO(ERROR) << "Workbench parser: something wrong happens while read EBLOCK.";
		}
	}

	for (size_t i = 0; i < _CMBlocks.size(); i++) {
		if (_CMBlocks[i].entity == CMBlock::Entity::NODE) {
			dMesh.nregions.push_back(MeshNRegion());
			dMesh.nregions.back().name = _CMBlocks[i].name;
			if (!_CMBlocks[i].readData(dMesh.nregions.back().nodes)) {
				ESINFO(ERROR) << "Workbench parser: something wrong happens while read CMBLOCK.";
			}
		} else {
			dMesh.eregions.push_back(MeshERegion());
			dMesh.eregions.back().name = _CMBlocks[i].name;
			if (!_CMBlocks[i].readData(dMesh.eregions.back().elements)) {
				ESINFO(ERROR) << "Workbench parser: something wrong happens while read CMBLOCK.";
			}
		}
	}
}
