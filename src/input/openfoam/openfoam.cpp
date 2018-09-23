
#include "openfoam.h"

#include "parser/points.h"
#include "parser/faces.h"
#include "parser/boundary.h"

#include "../plaindata.h"
#include "../../basis/containers/tarray.h"
#include "../../basis/logging/logging.h"
#include "../../basis/logging/timeeval.h"
#include "../../basis/utilities/communication.h"
#include "../../config/ecf/input/input.h"

#include "../../mesh/elements/element.h"
#include "../randominput.h"
#include "../sequentialinput.h"

#include <numeric>

using namespace espreso;

void OpenFOAMLoader::load(const InputConfiguration &configuration, Mesh &mesh)
{
	OpenFOAMLoader(configuration, mesh);
}

OpenFOAMLoader::OpenFOAMLoader(const InputConfiguration &configuration, Mesh &mesh)
: _configuration(configuration)
{
	TimeEval timing("Parsing OpenFOAM data");
	timing.totalTime.startWithBarrier();
	ESINFO(OVERVIEW) << "Load OpenFOAM data from directory '" << _configuration.path << "'.";

	TimeEvent tread("read data from file"); tread.start();
	readData();
	tread.end(); timing.addEvent(tread);
	ESINFO(PROGRESS2) << "OpenFOAM:: data copied from file.";

	PlainOpenFOAMData meshData;
	TimeEvent tparse("parsing data"); tparse.start();
	parseData(meshData);
	tparse.end(); timing.addEvent(tparse);
	ESINFO(PROGRESS2) << "OpenFOAM:: data parsed.";

	TimeEvent tbuild("building of elements"); tbuild.start();
	buildElements(meshData);
	tbuild.end(); timing.addEvent(tbuild);
	ESINFO(PROGRESS2) << "OpenFOAM:: elements builded.";

	TimeEvent tfaces("building of faces"); tfaces.start();
	buildFaces(meshData);
	tfaces.end(); timing.addEvent(tfaces);
	ESINFO(PROGRESS2) << "OpenFOAM:: faces builded.";

	timing.totalTime.endWithBarrier();
	timing.printStatsMPI();

	if (environment->MPIsize > 1) {
		RandomInput::buildMesh(meshData, mesh);
	} else {
		SequentialInput::buildMesh(meshData, mesh);
	}
}

void OpenFOAMLoader::readData()
{
	auto read = [&] (const std::string &file, ParallelFile &pfile) {
		if (!MPILoader::read(_configuration.path + "/constant/polyMesh/" + file, pfile, 80)) {
			ESINFO(ERROR) << "MPI cannot load file '" << _configuration.path + "/constant/polyMesh/" + file << "'";
		}
		MPILoader::align(pfile, 0);
	};

	read("points", _points);
	read("faces", _faces);
	read("neighbour", _neighbour);
	read("owner", _owner);
	read("boundary", _boundary);
}

void OpenFOAMLoader::parseData(PlainOpenFOAMData &mesh)
{
	if (!OpenFOAMPoints(_points.begin, _points.end).readData(mesh.nIDs, mesh.coordinates, _configuration.scale_factor)) {
		ESINFO(ERROR) << "OpenFOAM loader: cannot parse points.";
	}
	if (!OpenFOAMFaces(_faces.begin, _faces.end).readFaces(mesh)) {
		ESINFO(ERROR) << "OpenFOAM loader: cannot parse faces.";
	}
	if (!OpenFOAMFaces(_owner.begin, _owner.end).readParents(mesh.owner)) {
		ESINFO(ERROR) << "OpenFOAM loader: cannot parse owner.";
	}
	if (!OpenFOAMFaces(_neighbour.begin, _neighbour.end).readParents(mesh.neighbour)) {
		ESINFO(ERROR) << "OpenFOAM loader: cannot parse owner.";
	}

	if (!OpenFOAMBoundary(_boundary.begin, _boundary.end).readData(_boundaryData)) {
		ESINFO(ERROR) << "OpenFOAM loader: cannot parse boundary.";
	}
}

void OpenFOAMLoader::buildElements(PlainOpenFOAMData &mesh)
{
	size_t threads = environment->OMP_NUM_THREADS;

	auto sortIDs = [&] (std::vector<eslocal> &permutation, const std::vector<eslocal> &data) {
		permutation.resize(data.size());
		std::iota(permutation.begin(), permutation.end(), 0);

		std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) {
			if (data[i] == data[j]) {
				return i < j;
			}
			return data[i] < data[j];
		});
	};

	std::vector<eslocal> owner, neighbour;
	sortIDs(owner, mesh.owner);
	sortIDs(neighbour, mesh.neighbour);

	eslocal nelements = mesh.owner[owner.back()] + 1;
	std::vector<size_t> tdistribution = tarray<eslocal>::distribute(threads, nelements);

	std::vector<std::vector<eslocal> > esize(threads), enodes(threads);
	std::vector<std::vector<int> > etype(threads);

	_fdist.clear();
	_fdist.reserve(nelements + 1);
	_fdist.push_back(0);
	for (size_t f = 0; f < mesh.fsize.size(); f++) {
		_fdist.push_back(_fdist.back() + mesh.fsize[f]);
	}

	auto getThreadBegin = [] (const std::vector<eslocal> &data, const std::vector<eslocal> &perm, eslocal eindex) {
		return std::lower_bound(perm.begin(), perm.end(), eindex, [&] (eslocal i, eslocal eindex) {
			return data[perm[i]] < eindex;
		}) - perm.begin();
	};

	auto addFaces = [&] (const std::vector<eslocal> &data, const std::vector<eslocal> &perm, size_t &triangles, size_t &squares, size_t &index, eslocal element) {
		while (index < perm.size() && data[perm[index]] == element) {
			switch (mesh.fsize[perm[index++]]) {
			case 3: ++triangles; break;
			case 4:   ++squares; break;
			}
		}
	};

	auto getFace = [&] (const std::vector<eslocal> &data, const std::vector<eslocal> &perm, eslocal index, eslocal element, eslocal n1, eslocal n2) {
		while (index < perm.size() && data[perm[index]] == element) {
			for (eslocal f = 0; f < mesh.fsize[perm[index]]; f++) {
				if (n1 == mesh.fnodes[_fdist[perm[index]] + f] && n2 == mesh.fnodes[_fdist[perm[index]] + (f + 1) % mesh.fsize[perm[index]]]) {
					return std::pair<eslocal, eslocal>(index, f);
				}
			}
			++index;
		}
		return std::pair<eslocal, eslocal>(index, -1);
	};

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<eslocal> tsize, tnodes;
		std::vector<int> ttype;

		size_t oindex = getThreadBegin(mesh.owner, owner, tdistribution[t]);
		size_t nindex = getThreadBegin(mesh.neighbour, neighbour, tdistribution[t]);
		size_t triangles, squares;

		for (size_t e = tdistribution[t]; e < tdistribution[t + 1]; e++) {
			size_t obegin = oindex, nbegin = nindex, ebegin;
			std::pair<eslocal, eslocal> index;
			triangles = squares = 0;
			addFaces(mesh.owner, owner, triangles, squares, oindex, e);
			addFaces(mesh.neighbour, neighbour, triangles, squares, nindex, e);

			if (squares == 6 && triangles == 0) {
				ebegin = tnodes.size();
				ttype.push_back((int)Element::CODE::HEXA8);
				tsize.push_back(8);
				tnodes.insert(tnodes.end(), 8, -1);
				if (obegin < oindex) { // there is at least one owner
					tnodes[ebegin + 0] = mesh.fnodes[_fdist[owner[obegin]] + 0];
					tnodes[ebegin + 1] = mesh.fnodes[_fdist[owner[obegin]] + 1];
					tnodes[ebegin + 4] = mesh.fnodes[_fdist[owner[obegin]] + 3];
					tnodes[ebegin + 5] = mesh.fnodes[_fdist[owner[obegin]] + 2];
					++obegin;
				} else {
					tnodes[ebegin + 0] = mesh.fnodes[_fdist[neighbour[nbegin]] + 0];
					tnodes[ebegin + 1] = mesh.fnodes[_fdist[neighbour[nbegin]] + 3];
					tnodes[ebegin + 4] = mesh.fnodes[_fdist[neighbour[nbegin]] + 2];
					tnodes[ebegin + 5] = mesh.fnodes[_fdist[neighbour[nbegin]] + 1];
					++nbegin;
				}

				index = getFace(mesh.owner, owner, obegin, e, tnodes[ebegin + 1], tnodes[ebegin]);
				if (index.first < oindex) {
					tnodes[ebegin + 2] = mesh.fnodes[_fdist[owner[index.first]] + (index.second + 3) % mesh.fsize[index.first]];
					tnodes[ebegin + 3] = mesh.fnodes[_fdist[owner[index.first]] + (index.second + 2) % mesh.fsize[index.first]];
				} else {
					index = getFace(mesh.neighbour, neighbour, nbegin, e, tnodes[ebegin], tnodes[ebegin + 1]);
					tnodes[ebegin + 2] = mesh.fnodes[_fdist[neighbour[index.first]] + (index.second + 2) % mesh.fsize[index.first]];
					tnodes[ebegin + 3] = mesh.fnodes[_fdist[neighbour[index.first]] + (index.second + 3) % mesh.fsize[index.first]];
				}
				index = getFace(mesh.owner, owner, obegin, e, tnodes[ebegin + 4], tnodes[ebegin + 5]);
				if (index.first < oindex) {
					tnodes[ebegin + 6] = mesh.fnodes[_fdist[owner[index.first]] + (index.second + 2) % mesh.fsize[index.first]];
					tnodes[ebegin + 7] = mesh.fnodes[_fdist[owner[index.first]] + (index.second + 3) % mesh.fsize[index.first]];
				} else {
					index = getFace(mesh.neighbour, neighbour, nbegin, e, tnodes[ebegin + 5], tnodes[ebegin + 4]);
					tnodes[ebegin + 6] = mesh.fnodes[_fdist[neighbour[index.first]] + (index.second + 3) % mesh.fsize[index.first]];
					tnodes[ebegin + 7] = mesh.fnodes[_fdist[neighbour[index.first]] + (index.second + 2) % mesh.fsize[index.first]];
				}
				continue;
			}
			ESINFO(ERROR) << "OpenFOAM parser: an unknown element type with" << triangles << " triangles and " << squares << "squares.";
		}

		esize[t].swap(tsize);
		enodes[t].swap(tnodes);
		etype[t].swap(ttype);
	}

	for (size_t t = 0; t < threads; t++) {
		mesh.esize.insert(mesh.esize.end(), esize[t].begin(), esize[t].end());
		mesh.enodes.insert(mesh.enodes.end(), enodes[t].begin(), enodes[t].end());
		mesh.etype.insert(mesh.etype.end(), etype[t].begin(), etype[t].end());
	}
	mesh.body.resize(nelements, 0);
	mesh.material.resize(nelements, 0);
	mesh.eIDs.resize(nelements);
	std::iota(mesh.eIDs.begin(), mesh.eIDs.end(), 0);
}

void OpenFOAMLoader::buildFaces(PlainOpenFOAMData &mesh)
{
//	if (!Communication::broadcastUnknownSize(_boundaryData)) {
//		ESINFO(ERROR) << "ESPRESO internal error: exchange OpenFOAM boundary data.";
//	}

	for (size_t i = 0; i < _boundaryData.size(); i++) {
		mesh.esize.insert(mesh.esize.end(),
				mesh.fsize.begin() + _boundaryData[i].startFace,
				mesh.fsize.begin() + _boundaryData[i].startFace + _boundaryData[i].nFaces);
		mesh.enodes.insert(mesh.enodes.end(),
				mesh.fnodes.begin() + _fdist[_boundaryData[i].startFace],
				mesh.fnodes.begin() + _fdist[_boundaryData[i].startFace + _boundaryData[i].nFaces]);
		auto &indices = mesh.eregions[_boundaryData[i].name];
		indices.resize(_boundaryData[i].nFaces);
		std::iota(indices.begin(), indices.end(), mesh.esize.size() - _boundaryData[i].nFaces);
	}

	mesh.eIDs.reserve(mesh.esize.size());
	mesh.etype.reserve(mesh.esize.size());
	for (size_t e = mesh.eIDs.size(); e < mesh.esize.size(); e++) {
		mesh.eIDs.push_back(e);
		if (mesh.esize[e] == 3) {
			mesh.etype.push_back((int)Element::CODE::TRIANGLE3);
		}
		if (mesh.esize[e] == 4) {
			mesh.etype.push_back((int)Element::CODE::SQUARE4);
		}
	}
	mesh.body.resize(mesh.esize.size(), 0);
	mesh.material.resize(mesh.esize.size(), 0);
}
