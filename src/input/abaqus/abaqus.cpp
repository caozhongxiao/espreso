
#include "abaqus.h"

#include "parser/nblock.h"
#include "parser/eblock.h"
#include "parser/blockend.h"
#include "parser/nsel.h"
#include "parser/eset.h"
#include "parser/ssection.h"
#include "parser/material.h"
#include "parser/elemat.h"
#include "parser/nset.h"

#include "../../basis/containers/tarray.h"
#include "../../basis/logging/logging.h"
#include "../../basis/logging/timeeval.h"
#include "../../basis/utilities/communication.h"
#include "../../basis/utilities/utils.h"
#include "../../config/ecf/environment.h"
#include "../../config/ecf/input/input.h"

#include "../randominput.h"
#include "../sequentialinput.h"


#include "parser/parser.h"
// msg remove this

using namespace espreso;

void AbaqusLoader::load(const InputConfiguration &configuration, Mesh &mesh)
{
	AbaqusLoader(configuration, mesh);
}

AbaqusLoader::AbaqusLoader(const InputConfiguration &configuration, Mesh &mesh)
: _configuration(configuration)
{
	TimeEval timing("Parsing Workbench data");
	timing.totalTime.startWithBarrier();
	ESINFO(OVERVIEW) << "Load Abaqus data from '" << _configuration.path << "'.";

	TimeEvent tread("read data from file"); tread.start();
	readData();
	tread.end(); timing.addEvent(tread);
	ESINFO(PROGRESS2) << "Abaqus:: data copied from file.";

	TimeEvent tprepare("prepare data for parsing"); tprepare.start();
	prepareData();
	tprepare.end(); timing.addEvent(tprepare);
	ESINFO(PROGRESS2) << "Abaqus:: data prepared for parsing.";


	PlainAbaqusData meshData;
	TimeEvent tparse("parsing data"); tparse.start();
	parseData(meshData);
	tparse.end(); timing.addEvent(tparse);
	ESINFO(PROGRESS2) << "Abaqus:: data parsed.";
}

void AbaqusLoader::readData()
{
	TimeEval timing("Read data from file");
	timing.totalTime.startWithBarrier();

	MPISubset loaders(_configuration, MPITools::procs());

	TimeEvent e1("FILE OPEN");
	e1.start();

	MPI_File MPIFile;
	if (loaders.within.rank == 0) {
		if (!MPILoader::open(loaders.across, MPIFile, _configuration.path)) {
			ESINFO(ERROR) << "MPI cannot load file '" << _configuration.path << "'";
		}
	}

	e1.end();
	timing.addEvent(e1);

	TimeEvent e2("FILE READ");
	e2.start();

	if (loaders.within.rank == 0) {
		MPILoader::read(loaders.across, MPIFile, _pfile, MAX_LINE_STEP * MAX_LINE_SIZE);
	}

	e2.end();
	timing.addEvent(e2);

	TimeEvent e3("FILE SCATTER");
	e3.start();

	MPILoader::scatter(loaders.within, _pfile, MAX_LINE_STEP * MAX_LINE_SIZE);
	MPILoader::align(MPITools::procs(), _pfile, MAX_LINE_STEP);

	AbaqusParser::offset = _pfile.offsets[environment->MPIrank];
	AbaqusParser::begin = _pfile.begin;
	AbaqusParser::end = _pfile.end;

	e3.end();
	timing.addEvent(e3);

	timing.totalTime.endWithBarrier();
	timing.printStatsMPI();
}



void AbaqusLoader::prepareData()
{
	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<size_t> tdistribution = tarray<char>::distribute(threads, _pfile.end - _pfile.begin);
	std::vector<std::vector<NList> > tNLists(threads);
	std::vector<std::vector<EList> > tELists(threads);
	std::vector<std::vector<BlockFinish> > tBlockFinishs(threads);
	std::vector<std::vector<Eset> > tEsets(threads);
	std::vector<std::vector<Nset> > tNsets(threads);
	std::vector<std::vector<SSection> > tSSections(threads);
	std::vector<std::vector<Materials> > tMaterials(threads);

    #pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		const char* tbegin = _pfile.begin + tdistribution[t];
		if (t && *(tbegin - 1) != '\n') {
			while (tbegin < _pfile.end && *tbegin++ != '\n'); // start at new line
		}

		while (tbegin < _pfile.begin + tdistribution[t + 1]) {
			

			if (memcmp(tbegin, BlockFinish::asterik, BlockFinish::nSize) == 0  &&
					memcmp(tbegin , BlockFinish::comment, BlockFinish::ncSize) != 0) {
							tBlockFinishs[t].push_back(BlockFinish().parse(tbegin ));
						}

			if (memcmp(tbegin, NList::upper, NList::size) == 0) {
				tNLists[t].push_back(NList().parse(tbegin));
			}
			if (memcmp(tbegin, NList::lower, NList::size) == 0) {
				tNLists[t].push_back(NList().parse(tbegin));
			}
			if (memcmp(tbegin, NList::sentence, NList::size) == 0) {
				tNLists[t].push_back(NList().parse(tbegin));
	 		}
			


			if (memcmp(tbegin, Eset::upper, Eset::size) == 0) {
										tEsets[t].push_back(Eset().parse(tbegin));
			}
			if (memcmp(tbegin, Eset::lower, Eset::size) == 0) {
										tEsets[t].push_back(Eset().parse(tbegin));
			}
			if (memcmp(tbegin, Eset::sentence, Eset::size) == 0) {
							tEsets[t].push_back(Eset().parse(tbegin));
			}



		    if (memcmp(tbegin, EList::upper, EList::size) == 0) {
				tELists[t].push_back(EList().parse(tbegin));
			}
			if (memcmp(tbegin, EList::lower, EList::size) == 0) {
				tELists[t].push_back(EList().parse(tbegin));
			}
			if (memcmp(tbegin, EList::sentence, EList::size) == 0) {
				tELists[t].push_back(EList().parse(tbegin));
			}



		    if (memcmp(tbegin, SSection::upper, SSection::size) == 0) {
		    	tSSections[t].push_back(SSection().parse(tbegin));
			}
			if (memcmp(tbegin, SSection::lower, SSection::size) == 0) {
				tSSections[t].push_back(SSection().parse(tbegin));
			}
			if (memcmp(tbegin, SSection::sentence, SSection::size) == 0) {
				tSSections[t].push_back(SSection().parse(tbegin));
			}


		    if (memcmp(tbegin, Materials::upper, Materials::size) == 0) {
		    	tMaterials[t].push_back(Materials().parse(tbegin));
			}
			if (memcmp(tbegin, Materials::lower, Materials::size) == 0) {
				tMaterials[t].push_back(Materials().parse(tbegin));
			}
			if (memcmp(tbegin, Materials::sentence, Materials::size) == 0) {
				tMaterials[t].push_back(Materials().parse(tbegin));
			}



		    if (memcmp(tbegin, Nset::upper, Nset::size) == 0) {
		    	tNsets[t].push_back(Nset().parse(tbegin));
			}
			if (memcmp(tbegin, Nset::lower, Nset::size) == 0) {
				tNsets[t].push_back(Nset().parse(tbegin));
			}
			if (memcmp(tbegin, Nset::sentence, Nset::size) == 0) {
				tNsets[t].push_back(Nset().parse(tbegin));
			}

			while (tbegin < _pfile.end && *tbegin++ != '\n');
			
		}

	}
	for (size_t t = 0; t < threads; t++) {
		_NLists.insert(_NLists.end(), tNLists[t].begin(), tNLists[t].end());
		_blockFinishs.insert(_blockFinishs.end(), tBlockFinishs[t].begin(), tBlockFinishs[t].end());
		_ELists.insert(_ELists.end(), tELists[t].begin(), tELists[t].end());
		_Esets.insert(_Esets.end(), tEsets[t].begin(), tEsets[t].end());
		_SSections.insert(_SSections.end(), tSSections[t].begin(), tSSections[t].end());
		_Materials.insert(_Materials.end(), tMaterials[t].begin() ,tMaterials[t].end());
		_Nsets.insert(_Nsets.end(), tNsets[t].begin() ,tNsets[t].end());
	    }

	if (!Communication::allGatherUnknownSize(_blockFinishs)) {
						ESINFO(ERROR) << "ESPRESO internal error: exchange Workbench block ends.";
					}

	for (size_t i = 0; i < _NLists.size(); i++) {
		_NLists[i].fillDistribution(_blockFinishs, _pfile.offsets);
	    }

	for (size_t i = 0; i < _ELists.size(); i++) {
			_ELists[i].fillDistribution(_blockFinishs, _pfile.offsets);
		}

	for (size_t i = 0; i < _Esets.size(); i++) {
				_Esets[i].fillDistribution(_blockFinishs, _pfile.offsets);
		}

	for (size_t i = 0; i < _Nsets.size(); i++) {
				_Nsets[i].fillDistribution(_blockFinishs, _pfile.offsets);
		}

	if (!Communication::allGatherUnknownSize(_NLists)) {
						ESINFO(ERROR) << "ESPRESO internal error: exchange Workbench NBblocks.";
		}

	if (!Communication::allGatherUnknownSize(_ELists)) {
			ESINFO(ERROR) << "ESPRESO internal error: exchange Workbench EBlocks.";
		}

	if (!Communication::allGatherUnknownSize(_Esets)) {
				ESINFO(ERROR) << "ESPRESO internal error: exchange Workbench EBlocks.";
		}

	if (!Communication::allGatherUnknownSize(_Nsets)) {
				ESINFO(ERROR) << "ESPRESO internal error: exchange Workbench EBlocks.";
		}

	// fix distribution if EBlocks are across more processes and elements data have more lines
	for (size_t i = 0; i < _ELists.size(); i++) {
		_ELists[i].fixOffsets(_pfile.offsets);
		if (_pfile.begin != AbaqusParser::begin) {
			_pfile.offsets[environment->MPIrank] += _pfile.begin - AbaqusParser::begin;
		}
		if (_pfile.end != AbaqusParser::end) {
			_pfile.offsets[environment->MPIrank + 1] += _pfile.end - AbaqusParser::end;
		}
		_pfile.begin = AbaqusParser::begin;
		_pfile.end = AbaqusParser::end;
	}

}


void AbaqusLoader::parseData(PlainAbaqusData &meshData)
{

	//_Elemats[0].create_dict(_SSections,_Materials);

	for (size_t i = 0; i < 1;i++){   //_NLists.size(); i++) {
			if (!_NLists[i].readData(meshData.nIDs, meshData.coordinates, _configuration.scale_factor)) {
				ESINFO(ERROR) << "Abaqus parser: something wrong happens while read NList.";
			}
		}

	for (size_t i = 0; i < _ELists.size(); i++) {
			if (!_ELists[i].readData( meshData)) {
				ESINFO(ERROR) << "Workbench parser: something wrong happens while read EBLOCK.";
			}
		}

	for (size_t i = 0; i < _Esets.size(); i++) {
			if (!_Esets[i].readData( meshData.eregions[_Esets[i].NAME])) {
				ESINFO(ERROR) << "Workbench parser: something wrong happens while read EBLOCK.";
			}
		}

	for (size_t i = 0; i < _Nsets.size(); i++) {
			if (!_Nsets[i].readData( meshData.nregions[_Nsets[i].NAME])) {
				ESINFO(ERROR) << "Workbench parser: something wrong happens while read EBLOCK.";
			}
		}

	for (size_t i = 0; i < 1; i++) {

					ESINFO(ERROR) << "Workbench parser: something wrong happens while read EBLOCK.";

			}

}

