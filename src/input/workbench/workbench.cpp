
#include "workbench.h"

#include "parser/nblock.h"
#include "parser/eblock.h"
#include "parser/cmblock.h"
#include "parser/et.h"
#include "parser/esel.h"
#include "parser/nsel.h"
#include "parser/cm.h"
#include "parser/blockend.h"

#include "../../basis/containers/tarray.h"
#include "../../basis/logging/logging.h"
#include "../../basis/logging/timeeval.h"
#include "../../basis/utilities/communication.h"
#include "../../basis/utilities/utils.h"
#include "../../config/ecf/input/input.h"

#include "../randominput.h"
#include "../sequentialinput.h"

using namespace espreso;

void WorkbenchLoader::load(const InputConfiguration &configuration, Mesh &mesh)
{
	WorkbenchLoader(configuration, mesh);
}

WorkbenchLoader::WorkbenchLoader(const InputConfiguration &configuration, Mesh &mesh)
: _configuration(configuration)
{
	TimeEval timing("Parsing Workbench data");
	timing.totalTime.startWithBarrier();
	ESINFO(OVERVIEW) << "Load ANSYS Workbench data from '" << _configuration.path << "'.";

	TimeEvent tread("read data from file"); tread.start();
	readData();
	tread.end(); timing.addEvent(tread);
	ESINFO(PROGRESS2) << "Workbench:: data copied from file.";

	TimeEvent tprepare("prepare data for parsing"); tprepare.start();
	prepareData();
	tprepare.end(); timing.addEvent(tprepare);
	ESINFO(PROGRESS2) << "Workbench:: data prepared for parsing.";

	PlainWorkbenchData meshData;
	TimeEvent tparse("parsing data"); tparse.start();
	parseData(meshData);
	tparse.end(); timing.addEvent(tparse);
	ESINFO(PROGRESS2) << "Workbench:: data parsed.";

	timing.totalTime.endWithBarrier();
	timing.printStatsMPI();

	if (!_configuration.keep_material_sets) {
		std::fill(meshData.material.begin(), meshData.material.end(), 0);
	}

	if (environment->MPIsize > 1) {
		RandomInput::buildMesh(meshData, mesh);
	} else {
		SequentialInput::buildMesh(meshData, mesh);
	}
}

void WorkbenchLoader::readData()
{
	MPITools::MPICommunicator comm;
	Communication::createCommunicator(_configuration, comm);

	TimeEval timing("Read data from file");
	timing.totalTime.startWithBarrier();

	TimeEvent e1("FILE OPEN");
	e1.start();

	MPI_File MPIFile;
	if (comm.within.rank == 0) {
		if (!MPILoader::open(comm.across, MPIFile, _configuration.path)) {
			ESINFO(ERROR) << "MPI cannot load file '" << _configuration.path << "'";
		}
	}

	e1.end();
	timing.addEvent(e1);

	TimeEvent e2("FILE READ");
	e2.start();

	if (comm.within.rank == 0) {
		MPILoader::read(comm.across, MPIFile, _pfile, MAX_LINE_STEP * MAX_LINE_SIZE);
	}

	e2.end();
	timing.addEvent(e2);

	TimeEvent e3("FILE SCATTER");
	e3.start();

	MPILoader::scatter(comm.within, _pfile, MAX_LINE_STEP * MAX_LINE_SIZE);
	MPILoader::align(MPITools::procs(), _pfile, MAX_LINE_STEP);

	WorkbenchParser::offset = _pfile.offsets[environment->MPIrank];
	WorkbenchParser::begin = _pfile.begin;
	WorkbenchParser::end = _pfile.end;

	e3.end();
	timing.addEvent(e3);

	timing.totalTime.endWithBarrier();
	timing.printStatsMPI();
}

void WorkbenchLoader::prepareData()
{
	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<size_t> tdistribution = tarray<char>::distribute(threads, _pfile.end - _pfile.begin);
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
		const char* tbegin = _pfile.begin + tdistribution[t];
		if (t && *(tbegin - 1) != '\n') {
			while (tbegin < _pfile.end && *tbegin++ != '\n'); // start at new line
		}

		while (tbegin < _pfile.begin + tdistribution[t + 1]) {
			if (tbegin != _pfile.begin + tdistribution[t] && memcmp(tbegin - BlockEnd::unixSize, BlockEnd::unixEnd, BlockEnd::unixSize) == 0) {
				tBlockEnds[t].push_back(BlockEnd().parse(tbegin - BlockEnd::unixSize));
			}
			if (tbegin != _pfile.begin + tdistribution[t] && memcmp(tbegin - BlockEnd::winSize, BlockEnd::winEnd, BlockEnd::winSize) == 0) {
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
			while (tbegin < _pfile.end && *tbegin++ != '\n');
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
		_NBlocks[i].fillDistribution(_blockEnds, _pfile.offsets);
	}
	for (size_t i = 0; i < _EBlocks.size(); i++) {
		_EBlocks[i].fillDistribution(_blockEnds, _pfile.offsets);
	}
	for (size_t i = 0; i < _CMBlocks.size(); i++) {
		_CMBlocks[i].fillDistribution(_blockEnds, _pfile.offsets);
	}
	for (size_t i = 0; i < _ET.size(); i++) {
		_ET[i].fillDistribution(_blockEnds, _pfile.offsets);
	}
	for (size_t i = 0; i < _ESel.size(); i++) {
		_ESel[i].fillDistribution(_blockEnds, _pfile.offsets);
	}
	for (size_t i = 0; i < _NSel.size(); i++) {
		_NSel[i].fillDistribution(_blockEnds, _pfile.offsets);
	}
	for (size_t i = 0; i < _CM.size(); i++) {
		_CM[i].fillDistribution(_blockEnds, _pfile.offsets);
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
		_EBlocks[i].fixOffsets(_pfile.offsets);
		if (_pfile.begin != WorkbenchParser::begin) {
			_pfile.offsets[environment->MPIrank] += _pfile.begin - WorkbenchParser::begin;
		}
		if (_pfile.end != WorkbenchParser::end) {
			_pfile.offsets[environment->MPIrank + 1] += _pfile.end - WorkbenchParser::end;
		}
		_pfile.begin = WorkbenchParser::begin;
		_pfile.end = WorkbenchParser::end;
	}
}

void WorkbenchLoader::parseData(PlainWorkbenchData &meshData)
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
		if (!_NBlocks[i].readData(meshData.nIDs, meshData.coordinates, _configuration.scale_factor)) {
			ESINFO(ERROR) << "Workbench parser: something wrong happens while read NBLOCK.";
		}
	}

	for (size_t i = 0; i < _EBlocks.size(); i++) {
		if (!_EBlocks[i].readData(_ET, meshData)) {
			ESINFO(ERROR) << "Workbench parser: something wrong happens while read EBLOCK.";
		}
	}

	for (size_t i = 0; i < _CMBlocks.size(); i++) {
		switch (_CMBlocks[i].entity) {
		case CMBlock::Entity::NODE: {
			if (!_CMBlocks[i].readData(meshData.nregions[_CMBlocks[i].name])) {
				ESINFO(ERROR) << "Workbench parser: something wrong happens while read CMBLOCK.";
			}
		} break;
		case CMBlock::Entity::ELEMENT: {
			if (!_CMBlocks[i].readData(meshData.eregions[_CMBlocks[i].name])) {
				ESINFO(ERROR) << "Workbench parser: something wrong happens while read CMBLOCK.";
			}
		} break;
		default:
			ESINFO(ERROR) << "ESPRESO Workbench parser: unknown CMBLOCK type.";
		}
	}

	for (size_t i = 0; i < _CM.size(); i++) {
		_CM[i].addRegion(meshData, _ESel, meshData.eregions, _NSel, meshData.nregions);
	}
}
