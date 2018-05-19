
#include "workbench.h"
#include "../../basis/containers/tarray.h"
#include "../../basis/logging/logging.h"
#include "../../basis/logging/timeeval.h"
#include "../../basis/utilities/communication.h"
#include "../../basis/utilities/utils.h"
#include "../../config/ecf/root.h"

#include "../../mesh/mesh.h"
#include "../randominput.h"

using namespace espreso;

void WorkbenchLoader::load(const ECFRoot &configuration, Mesh &mesh)
{
	WorkbenchLoader(configuration, mesh);
}

WorkbenchLoader::WorkbenchLoader(const ECFRoot &configuration, Mesh &mesh)
: _configuration(configuration), _mesh(mesh)
{
	TimeEval timing("Parsing Workbench data");
	timing.totalTime.startWithBarrier();
	ESINFO(OVERVIEW) << "Load ANSYS Workbench data from '" << configuration.workbench.path << "'.";

	TimeEvent tread("read data from file"); tread.start();
	readData();
	tread.end(); timing.addEvent(tread);
	ESINFO(PROGRESS2) << "Workbench:: data copied from file.";

	TimeEvent tprepare("prepare data for parsing"); tprepare.start();
	prepareData();
	tprepare.end(); timing.addEvent(tprepare);
	ESINFO(PROGRESS2) << "Workbench:: data prepared for parsing.";

	PlainMeshData meshData;
	TimeEvent tparse("parsing data"); tparse.start();
	parseData(meshData);
	tparse.end(); timing.addEvent(tparse);
	ESINFO(PROGRESS2) << "Workbench:: data parsed.";

	timing.totalTime.endWithBarrier();
	timing.printStatsMPI();

	RandomInput::buildMesh(configuration, meshData, mesh);
}

void WorkbenchLoader::readData()
{
	TimeEval timing("Read data from file");
	timing.totalTime.startWithBarrier();

	TimeEvent e1("prepare");
	e1.start();

	MPI_File MPIfile;

	if (MPI_File_open(environment->MPICommunicator, _configuration.workbench.path.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &MPIfile)) {
		ESINFO(ERROR) << "MPI cannot load file '" << _configuration.workbench.path << "'";
	}

	MPI_Offset size;
	MPI_File_get_size(MPIfile, &size);

	size_t block = 1;
	while (size / (1L << (block - 1)) > (1L << 31)) {
		++block;
	}
	block = 1L << block;

	std::vector<size_t> fdistribution = tarray<int>::distribute(environment->MPIsize, size / block + ((size % block) ? 1 : 0));
	std::vector<MPI_Aint> displacement = { (MPI_Aint)(block * fdistribution[environment->MPIrank]) };
	std::vector<int> length = { (int)(fdistribution[environment->MPIrank + 1] - fdistribution[environment->MPIrank]) };

	MPI_Datatype chunk;
	MPI_Datatype fDataDistribution;

	MPI_Type_contiguous(block, MPI_BYTE, &chunk);
	MPI_Type_commit(&chunk);
	MPI_Type_create_hindexed(1, length.data(), displacement.data(), chunk, &fDataDistribution);
	MPI_Type_commit(&fDataDistribution);

	_data.resize(block * length.front() + MAX_LINE_STEP * MAX_LINE_SIZE);

	MPI_File_set_view(MPIfile, 0, chunk, fDataDistribution, "native", MPI_INFO_NULL);

	e1.end();
	timing.addEvent(e1);

	TimeEvent e2("read file");
	e2.start();

	MPI_File_read_all(MPIfile, _data.data(), length.front(), chunk, MPI_STATUS_IGNORE);

	e2.end();
	timing.addEvent(e2);

	TimeEvent e3("post-process");
	e3.start();

	MPI_File_close(&MPIfile);

	MPI_Type_free(&chunk);

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
	if (environment->MPIrank + 1 == environment->MPIsize) {
		_end -= block * fdistribution.back() - size;
	}
	_begin = _current;
	_dataOffset = Communication::getDistribution<size_t>(_end - _current, MPITools::operations().sizeToOffsetsSize_t);

	WorkbenchParser::offset = _dataOffset[environment->MPIrank];
	WorkbenchParser::begin = _begin;
	WorkbenchParser::end = _end;

	e3.end();
	timing.addEvent(e3);

	timing.totalTime.endWithBarrier();
	timing.printStatsMPI();
}

void WorkbenchLoader::prepareData()
{
	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<size_t> tdistribution = tarray<char>::distribute(threads, _end - _begin);
	std::vector<std::vector<NBlock> > tNBlocks(threads);
	std::vector<std::vector<EBlock> > tEBlocks(threads);
	std::vector<std::vector<CMBlock> > tCMBlocks(threads);
	std::vector<std::vector<ET> > tET(threads);
	std::vector<std::vector<ESel> > tESel(threads);
	std::vector<std::vector<NSel> > tNSel(threads);
	std::vector<std::vector<CM> > tCM(threads);
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
			if (memcmp(tbegin, ET::upper, ET::size) == 0) {
				tET[t].push_back(ET().parse(tbegin));
			}
			if (memcmp(tbegin, ET::lower, ET::size) == 0) {
				tET[t].push_back(ET().parse(tbegin));
			}
			if (memcmp(tbegin, ESel::upper, ESel::size) == 0) {
				tESel[t].push_back(ESel().parse(tbegin));
			}
			if (memcmp(tbegin, ESel::lower, ESel::size) == 0) {
				tESel[t].push_back(ESel().parse(tbegin));
			}
			// TODO: implement
//			if (memcmp(tbegin, NSel::upper, NSel::size) == 0) {
//				tNSel[t].push_back(NSel().parse(tbegin));
//			}
//			if (memcmp(tbegin, NSel::lower, NSel::size) == 0) {
//				tNSel[t].push_back(NSel().parse(tbegin));
//			}
			if (memcmp(tbegin, CM::upper, CM::size) == 0) {
				tCM[t].push_back(CM().parse(tbegin));
			}
			if (memcmp(tbegin, CM::lower, CM::size) == 0) {
				tCM[t].push_back(CM().parse(tbegin));
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
		_ET.insert(_ET.end(), tET[t].begin(), tET[t].end());
		_ESel.insert(_ESel.end(), tESel[t].begin(), tESel[t].end());
		_NSel.insert(_NSel.end(), tNSel[t].begin(), tNSel[t].end());
		_CM.insert(_CM.end(), tCM[t].begin(), tCM[t].end());
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
	for (size_t i = 0; i < _ET.size(); i++) {
		_ET[i].fillDistribution(_blockEnds, _dataOffset);
	}
	for (size_t i = 0; i < _ESel.size(); i++) {
		_ESel[i].fillDistribution(_blockEnds, _dataOffset);
	}
	for (size_t i = 0; i < _NSel.size(); i++) {
		_NSel[i].fillDistribution(_blockEnds, _dataOffset);
	}
	for (size_t i = 0; i < _CM.size(); i++) {
		_CM[i].fillDistribution(_blockEnds, _dataOffset);
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

	if (!Communication::allGatherUnknownSize(_ET)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange Workbench ET.";
	}
	if (!Communication::allGatherUnknownSize(_ESel)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange Workbench ESel.";
	}
	if (!Communication::allGatherUnknownSize(_NSel)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange Workbench NSel.";
	}
	if (!Communication::allGatherUnknownSize(_CM)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange Workbench CM.";
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
//	for (size_t i = _EBlocks.size() - 1; i < _EBlocks.size(); i--) {
//		if (!_EBlocks[i].Solkey) {
//			_BBlocks.push_back(_EBlocks[i]);
//			_EBlocks.erase(_EBlocks.begin() + i);
//		}
//	}
}

void WorkbenchLoader::parseData(PlainMeshData &meshData)
{
	std::vector<ET> et;
	eslocal maxet = 0;
	for (size_t i = 0; i < _ET.size(); i++) {
		if (maxet < _ET[i].id) {
			maxet = _ET[i].id;
		}
	}
	et.resize(maxet + 1);
	for (size_t i = 0; i < _ET.size(); i++) {
		if (_ET[i].id >= 0) {
			et[_ET[i].id] = _ET[i];
		}
	}
	_ET.swap(et);

	for (size_t i = 0; i < _NBlocks.size(); i++) {
		if (!_NBlocks[i].readData(meshData.nIDs, meshData.coordinates, _configuration.workbench.scale_factor)) {
			ESINFO(ERROR) << "Workbench parser: something wrong happens while read NBLOCK.";
		}
	}

	for (size_t i = 0; i < _EBlocks.size(); i++) {
		if (!_EBlocks[i].readData(_ET, meshData.esize, meshData.enodes, meshData.edata)) {
			ESINFO(ERROR) << "Workbench parser: something wrong happens while read EBLOCK.";
		}
	}

//	for (size_t i = 0; i < _BBlocks.size(); i++) {
//		meshData.bregions.push_back(MeshBRegion());
//		meshData.bregions.back().name = "NAMELESS_SET_" + std::to_string(i + 1);
//		if (!_BBlocks[i].readData(_ET, meshData.bregions.back().esize, meshData.bregions.back().enodes, meshData.bregions.back().edata)) {
//			ESINFO(ERROR) << "Workbench parser: something wrong happens while read EBLOCK.";
//		}
//		meshData.bregions.back().min = meshData.bregions.back().edata.size() ? meshData.bregions.back().edata.front().id : 0;
//		meshData.bregions.back().max = meshData.bregions.back().edata.size() ? meshData.bregions.back().edata.back().id : 0;
//		MPI_Bcast(&meshData.bregions.back().min, sizeof(eslocal), MPI_BYTE, _BBlocks[i].fRank, environment->MPICommunicator);
//		MPI_Bcast(&meshData.bregions.back().max, sizeof(eslocal), MPI_BYTE, _BBlocks[i].lRank, environment->MPICommunicator);
//	}

	for (size_t i = 0; i < _CMBlocks.size(); i++) {
		switch (_CMBlocks[i].entity) {
		case CMBlock::Entity::NODE: {
			meshData.nregions.push_back(MeshNRegion());
			meshData.nregions.back().name = _CMBlocks[i].name;
			if (!_CMBlocks[i].readData(meshData.nregions.back().nodes)) {
				ESINFO(ERROR) << "Workbench parser: something wrong happens while read CMBLOCK.";
			}
		} break;
		case CMBlock::Entity::ELEMENT: {
			meshData.eregions.push_back(MeshERegion());
			meshData.eregions.back().name = _CMBlocks[i].name;
			if (!_CMBlocks[i].readData(meshData.eregions.back().elements)) {
				ESINFO(ERROR) << "Workbench parser: something wrong happens while read CMBLOCK.";
			}
			meshData.eregions.back().min = meshData.eregions.back().elements.size() ? meshData.eregions.back().elements.front() : 0;
			meshData.eregions.back().max = meshData.eregions.back().elements.size() ? meshData.eregions.back().elements.back() : 0;
			MPI_Bcast(&meshData.eregions.back().min, sizeof(eslocal), MPI_BYTE, _CMBlocks[i].fRank, environment->MPICommunicator);
			MPI_Bcast(&meshData.eregions.back().max, sizeof(eslocal), MPI_BYTE, _CMBlocks[i].lRank, environment->MPICommunicator);
		} break;
		default:
			ESINFO(ERROR) << "ESPRESO Workbench parser: unknown CMBLOCK type.";
		}
	}

	for (size_t i = 0; i < _CM.size(); i++) {
		_CM[i].addRegion(_ESel, meshData.edata, meshData.eregions, _NSel, meshData.nregions);
		meshData.eregions.back().min = meshData.eregions.back().elements.size() ? meshData.eregions.back().elements.front() : 0;
		meshData.eregions.back().max = meshData.eregions.back().elements.size() ? meshData.eregions.back().elements.back() : 0;
		MPI_Bcast(&meshData.eregions.back().min, sizeof(eslocal), MPI_BYTE, _CM[i].fRank, environment->MPICommunicator);
		MPI_Bcast(&meshData.eregions.back().max, sizeof(eslocal), MPI_BYTE, _CM[i].lRank, environment->MPICommunicator);
	}
}
