
#include "openfoam.h"

#include "parser/points.h"
#include "parser/faces.h"
#include "parser/boundary.h"

#include "../plaindata.h"
#include "../../basis/containers/tarray.h"
#include "../../basis/logging/logging.h"
#include "../../basis/logging/timeeval.h"
#include "../../basis/utilities/utils.h"
#include "../../basis/utilities/communication.h"
#include "../../config/ecf/environment.h"
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

	TimeEvent tcollect("collection of faces"); tcollect.start();
	collectFaces(meshData);
	tcollect.end(); timing.addEvent(tcollect);
	ESINFO(PROGRESS2) << "OpenFOAM:: faces collected.";

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
	ProcessesReduction onlyroot;
	onlyroot.granularity = ProcessesReduction::Granularity::PROCESSES;
	onlyroot.pattern = ProcessesReduction::Pattern::PREFIX;
	onlyroot.reduction_ratio = 1;

	MPISubset subset(_configuration, MPITools::procs()), singleton(onlyroot, MPITools::procs());

	auto distributedReader = [&] (const std::string &file, ParallelFile &pfile) {
		if (subset.within.rank == 0) {
			MPI_File MPIFile;
			if (!MPILoader::open(subset.across, MPIFile, _configuration.path + "/constant/polyMesh/" + file)) {
				ESINFO(ERROR) << "MPI cannot load file '" << _configuration.path + "/constant/polyMesh/" + file << "'";
			}
			MPILoader::read(subset.across, MPIFile, pfile, 80);
		}

		MPILoader::scatter(subset.within, pfile, 80);
		MPILoader::align(MPITools::procs(), pfile, 0);
	};

	auto singletonReader = [&] (const std::string &file, ParallelFile &pfile) {
		if (singleton.within.rank == 0) {
			MPI_File MPIFile;
			if (!MPILoader::open(singleton.across, MPIFile, _configuration.path + "/constant/polyMesh/" + file)) {
				ESINFO(ERROR) << "MPI cannot load file '" << _configuration.path + "/constant/polyMesh/" + file << "'";
			}
			MPILoader::read(singleton.across, MPIFile, pfile, 80);
		}
		MPILoader::bcast(singleton.within, pfile);
	};

	distributedReader("points", _points);
	distributedReader("faces", _faces);
	distributedReader("neighbour", _neighbour);
	distributedReader("owner", _owner);
	singletonReader("boundary", _boundary);
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

void OpenFOAMLoader::collectFaces(PlainOpenFOAMData &mesh)
{
	std::vector<size_t> ownersDist = Communication::getDistribution(mesh.owner.size());
	std::vector<size_t> neighborsDist = Communication::getDistribution(mesh.neighbour.size());
	std::vector<size_t> target = Communication::getDistribution(mesh.fsize.size());

	if (!Communication::balance(mesh.owner, ownersDist, target)) {
		ESINFO(ERROR) << "ESPRESO internal error: balance faces owners.";
	}

	for (size_t i = 0; i < target.size(); i++) {
		if (target[i] > neighborsDist.back()) {
			target[i] = neighborsDist.back();
		}
	}

	if (!Communication::balance(mesh.neighbour, neighborsDist, target)) {
		ESINFO(ERROR) << "ESPRESO internal error: balance faces owners.";
	}

	eslocal maxID = 0, gmaxID;
	if (mesh.owner.size()) {
		maxID = *std::max_element(mesh.owner.begin(), mesh.owner.end());
	}
	if (mesh.neighbour.size()) {
		maxID = std::max(*std::max_element(mesh.neighbour.begin(), mesh.neighbour.end()), maxID);
	}

	MPI_Allreduce(&maxID, &gmaxID, sizeof(eslocal), MPI_BYTE, MPITools::eslocalOperations().max, environment->MPICommunicator);

	_edist = tarray<eslocal>::distribute(environment->MPIsize, gmaxID + 1);
	size_t firstID = Communication::getDistribution(mesh.fsize.size())[environment->MPIrank];

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

	std::vector<eslocal> oPermutation, nPermutation;
	sortIDs(oPermutation, mesh.owner);
	sortIDs(nPermutation, mesh.neighbour);

	std::vector<eslocal> fdist({ 0 });
	fdist.reserve(mesh.fsize.size() + 1);
	for (size_t e = 0; e < mesh.fsize.size(); e++) {
		fdist.push_back(fdist.back() + mesh.fsize[e]);
	}

	std::vector<eslocal> sBuffer, rBuffer;
	// ID, size, owner / -1 * neighbor, nodes
	sBuffer.reserve(4 * environment->MPIsize + 3 * (mesh.owner.size() + mesh.neighbour.size()) + mesh.fnodes.size());

	size_t prevsize;
	auto obegin = oPermutation.begin();
	auto nbegin = nPermutation.begin();
	for (int r = 0; r < environment->MPIsize; r++) {
		prevsize = sBuffer.size();
		sBuffer.push_back(0); // total size
		sBuffer.push_back(r); // target
		sBuffer.push_back(0); // number of faces

		auto o = obegin;
		for ( ; o != oPermutation.end() && mesh.owner[*o] < _edist[r + 1]; ++o) {
			sBuffer.push_back(firstID + *o);
			sBuffer.push_back(mesh.owner[*o]);
			sBuffer.push_back(mesh.fsize[*o]);
			sBuffer.insert(sBuffer.end(), mesh.fnodes.begin() + fdist[*o], mesh.fnodes.begin() + fdist[*o + 1]);
		}
		auto n = nbegin;
		for ( ; n != nPermutation.end() && mesh.neighbour[*n] < _edist[r + 1]; ++n) {
			sBuffer.push_back(firstID + *n);
			sBuffer.push_back(-mesh.neighbour[*n] - 1);
			sBuffer.push_back(mesh.fsize[*n]);
			sBuffer.insert(sBuffer.end(), mesh.fnodes.begin() + fdist[*n], mesh.fnodes.begin() + fdist[*n + 1]);
		}
		sBuffer[prevsize + 2] = (o - obegin) + (n - nbegin);
		obegin = o;
		nbegin = n;

		sBuffer[prevsize] = sBuffer.size() - prevsize;
	}

	if (!Communication::allToAllWithDataSizeAndTarget(sBuffer, rBuffer)) {
		ESINFO(ERROR) << "ESPRESO internal error: distribute permuted elements.";
	}

	mesh.fsize.clear();
	mesh.fnodes.clear();
	mesh.owner.clear();
	mesh.neighbour.clear();

	std::vector<eslocal> fIDs, fsize, fnodes, owners;

	size_t offset = 0;
	for (int r = 0; r < environment->MPIsize; r++) {
		++offset;
		size_t size = rBuffer[++offset];
		++offset;

		for (size_t f = 0; f < size; ++f) {
			fIDs.push_back(rBuffer[offset++]);
			owners.push_back(rBuffer[offset++]);
			fsize.push_back(rBuffer[offset++]);
			fnodes.insert(fnodes.end(), rBuffer.begin() + offset, rBuffer.begin() + offset + fsize.back());
			offset += fsize.back();
		}
	}

	std::vector<eslocal> fpermutation(fIDs.size());
	std::iota(fpermutation.begin(), fpermutation.end(), 0);
	std::sort(fpermutation.begin(), fpermutation.end(), [&] (eslocal i, eslocal j) {
		return fIDs[i] < fIDs[j];
	});

	fdist.resize(1);
	fdist.reserve(fsize.size() + 1);
	for (size_t e = 0; e < fsize.size(); e++) {
		fdist.push_back(fdist.back() + fsize[e]);
	}

	if (fpermutation.size()) {
		size_t rest = 0;
		for (auto i = fpermutation.begin(); i != fpermutation.end() - 1; ++i) {
			if (fIDs[*i] == fIDs[*(i + 1)]) {
				mesh.fsize.push_back(fsize[*i]);
				mesh.fID.push_back(fIDs[*i]);
				if (owners[*i] >= 0) {
					mesh.owner.push_back(owners[*i]);
					mesh.neighbour.push_back(-owners[*(i + 1)] - 1);
					mesh.fnodes.insert(mesh.fnodes.end(), fnodes.begin() + fdist[*i], fnodes.begin() + fdist[*i + 1]);
				} else {
					mesh.owner.push_back(owners[*(i + 1)]);
					mesh.neighbour.push_back(-owners[*i] - 1);
					mesh.fnodes.insert(mesh.fnodes.end(), fnodes.begin() + fdist[*(i + 1)], fnodes.begin() + fdist[*(i + 1) + 1]);
				}
				++i;
			} else {
				fpermutation[rest++] = *i;
			}
		}
		if (fIDs[*(fpermutation.end() - 1)] != fIDs[*(fpermutation.end() - 2)]) {
			fpermutation[rest++] = fpermutation.back();
		}

		for (auto i = fpermutation.begin(); i != fpermutation.begin() + rest; ++i) {
			mesh.fsize.push_back(fsize[*i]);
			mesh.fID.push_back(fIDs[*i]);
			if (owners[*i] >= 0) {
				mesh.owner.push_back(owners[*i]);
				mesh.fnodes.insert(mesh.fnodes.end(), fnodes.begin() + fdist[*i], fnodes.begin() + fdist[*i + 1]);
			} else {
				mesh.owner.push_back(-owners[*i] - 1);
				mesh.fnodes.insert(mesh.fnodes.end(), fnodes.rbegin() + fdist.back() - fdist[*i + 1], fnodes.rbegin() + fdist.back() - fdist[*i]);
			}
		}
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

	std::vector<size_t> tdistribution = tarray<size_t>::distribute(threads, _edist[environment->MPIrank + 1] - _edist[environment->MPIrank]);

	std::vector<std::vector<eslocal> > esize(threads), enodes(threads);
	std::vector<std::vector<int> > etype(threads);

	_fdist.clear();
	_fdist.reserve(tdistribution.back() + 1);
	_fdist.push_back(0);
	for (size_t f = 0; f < mesh.fsize.size(); f++) {
		_fdist.push_back(_fdist.back() + mesh.fsize[f]);
	}

	auto getThreadBegin = [] (const std::vector<eslocal> &data, const std::vector<eslocal> &perm, eslocal eindex) {
		return std::lower_bound(perm.begin(), perm.end(), eindex, [&] (eslocal i, eslocal eindex) {
			return data[i] < eindex;
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

	auto getUnknown = [&] (eslocal *kbegin, eslocal *kend, eslocal *ubegin, eslocal *uend) {
		for (auto i = ubegin, j = kbegin; i != uend; ++i) {
			for (j = kbegin; j != kend; ++j) {
				if (*i == *j) {
					break;
				}
			}
			if (j == kend) {
				return *i;
			}
		}
	};

	auto findElementWithSize = [&] (const std::vector<eslocal> &perm, size_t &index, size_t &max, size_t size) {
		while (index < max) { // there is at least one owner
			if (mesh.fsize[perm[index]] == size) {
				break;
			} else {
				++index;
			}
		}
	};

	size_t eoffset = _edist[environment->MPIrank];

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<eslocal> tsize, tnodes;
		std::vector<int> ttype;

		size_t oindex = getThreadBegin(mesh.owner, owner, tdistribution[t] + eoffset);
		size_t nindex = getThreadBegin(mesh.neighbour, neighbour, tdistribution[t] + eoffset);

		for (size_t e = tdistribution[t] + eoffset; e < tdistribution[t + 1] + eoffset; e++) {
			size_t obegin = oindex, nbegin = nindex;
			std::pair<eslocal, eslocal> index;
			size_t triangles = 0, squares = 0;
			addFaces(mesh.owner, owner, triangles, squares, oindex, e);
			addFaces(mesh.neighbour, neighbour, triangles, squares, nindex, e);

			if (squares == 6 && triangles == 0) {
				size_t ebegin = tnodes.size();
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
					tnodes[ebegin + 2] = mesh.fnodes[_fdist[owner[index.first]] + (index.second + 3) % mesh.fsize[owner[index.first]]];
					tnodes[ebegin + 3] = mesh.fnodes[_fdist[owner[index.first]] + (index.second + 2) % mesh.fsize[owner[index.first]]];
				} else {
					index = getFace(mesh.neighbour, neighbour, nbegin, e, tnodes[ebegin], tnodes[ebegin + 1]);
					tnodes[ebegin + 2] = mesh.fnodes[_fdist[neighbour[index.first]] + (index.second + 2) % mesh.fsize[neighbour[index.first]]];
					tnodes[ebegin + 3] = mesh.fnodes[_fdist[neighbour[index.first]] + (index.second + 3) % mesh.fsize[neighbour[index.first]]];
				}
				index = getFace(mesh.owner, owner, obegin, e, tnodes[ebegin + 4], tnodes[ebegin + 5]);
				if (index.first < oindex) {
					tnodes[ebegin + 6] = mesh.fnodes[_fdist[owner[index.first]] + (index.second + 2) % mesh.fsize[owner[index.first]]];
					tnodes[ebegin + 7] = mesh.fnodes[_fdist[owner[index.first]] + (index.second + 3) % mesh.fsize[owner[index.first]]];
				} else {
					index = getFace(mesh.neighbour, neighbour, nbegin, e, tnodes[ebegin + 5], tnodes[ebegin + 4]);
					tnodes[ebegin + 6] = mesh.fnodes[_fdist[neighbour[index.first]] + (index.second + 3) % mesh.fsize[neighbour[index.first]]];
					tnodes[ebegin + 7] = mesh.fnodes[_fdist[neighbour[index.first]] + (index.second + 2) % mesh.fsize[neighbour[index.first]]];
				}
				continue;
			}

			if (squares == 0 && triangles == 4) {
				size_t ebegin = tnodes.size();
				ttype.push_back((int)Element::CODE::TETRA4);
				tsize.push_back(4);
				tnodes.insert(tnodes.end(), 4, -1);
				if (obegin < oindex) {
					tnodes[ebegin + 0] = mesh.fnodes[_fdist[owner[obegin]] + 0];
					tnodes[ebegin + 1] = mesh.fnodes[_fdist[owner[obegin]] + 2];
					tnodes[ebegin + 2] = mesh.fnodes[_fdist[owner[obegin]] + 1];
					++obegin;
				} else {
					tnodes[ebegin + 0] = mesh.fnodes[_fdist[neighbour[nbegin]] + 0];
					tnodes[ebegin + 1] = mesh.fnodes[_fdist[neighbour[nbegin]] + 1];
					tnodes[ebegin + 2] = mesh.fnodes[_fdist[neighbour[nbegin]] + 2];
					++nbegin;
				}

				if (obegin < oindex) {
					tnodes[ebegin + 3] = getUnknown(
							tnodes.data() + ebegin, tnodes.data() + ebegin + 3,
							mesh.fnodes.data() + _fdist[owner[obegin]], mesh.fnodes.data() + _fdist[owner[obegin] + 1]);
				} else {
					tnodes[ebegin + 3] = getUnknown(
							tnodes.data() + ebegin, tnodes.data() + ebegin + 3,
							mesh.fnodes.data() + _fdist[neighbour[nbegin]], mesh.fnodes.data() + _fdist[neighbour[nbegin] + 1]);
				}
				continue;
			}

			if (squares == 3 && triangles == 2) {
				size_t ebegin = tnodes.size();
				ttype.push_back((int)Element::CODE::PRISMA6);
				tsize.push_back(6);
				tnodes.insert(tnodes.end(), 6, -1);
				size_t otria = obegin, ntria = nbegin;
				findElementWithSize(owner, otria, oindex, 3);
				findElementWithSize(neighbour, ntria, nindex, 3);

				if (otria < oindex) {
					tnodes[ebegin + 0] = mesh.fnodes[_fdist[owner[otria]] + 0];
					tnodes[ebegin + 1] = mesh.fnodes[_fdist[owner[otria]] + 1];
					tnodes[ebegin + 2] = mesh.fnodes[_fdist[owner[otria]] + 2];
					findElementWithSize(owner, ++otria, oindex, 3);
				} else {
					tnodes[ebegin + 0] = mesh.fnodes[_fdist[neighbour[ntria]] + 0];
					tnodes[ebegin + 1] = mesh.fnodes[_fdist[neighbour[ntria]] + 2];
					tnodes[ebegin + 2] = mesh.fnodes[_fdist[neighbour[ntria]] + 1];
					findElementWithSize(neighbour, ++ntria, nindex, 3);
				}

				index = getFace(mesh.owner, owner, obegin, e, tnodes[ebegin + 1], tnodes[ebegin]);
				if (mesh.fsize[owner[index.first]] == 3) {
					index = getFace(mesh.owner, owner, ++obegin, e, tnodes[ebegin + 1], tnodes[ebegin]);
				}
				if (index.first < oindex) {
					tnodes[ebegin + 3] = mesh.fnodes[_fdist[owner[index.first]] + (index.second + 2) % mesh.fsize[owner[index.first]]];
					tnodes[ebegin + 4] = mesh.fnodes[_fdist[owner[index.first]] + (index.second + 3) % mesh.fsize[owner[index.first]]];
				} else {
					index = getFace(mesh.neighbour, neighbour, nbegin, e, tnodes[ebegin], tnodes[ebegin + 1]);
					if (mesh.fsize[neighbour[index.first]] == 3) {
						index = getFace(mesh.neighbour, neighbour, ++nbegin, e, tnodes[ebegin], tnodes[ebegin + 1]);
					}
					tnodes[ebegin + 3] = mesh.fnodes[_fdist[neighbour[index.first]] + (index.second + 3) % mesh.fsize[neighbour[index.first]]];
					tnodes[ebegin + 4] = mesh.fnodes[_fdist[neighbour[index.first]] + (index.second + 2) % mesh.fsize[neighbour[index.first]]];
				}

				if (otria < oindex) {
					tnodes[ebegin + 5] = getUnknown(
							tnodes.data() + ebegin, tnodes.data() + ebegin + 5,
							mesh.fnodes.data() + _fdist[owner[otria]], mesh.fnodes.data() + _fdist[owner[otria] + 1]);
				} else {
					tnodes[ebegin + 5] = getUnknown(
							tnodes.data() + ebegin, tnodes.data() + ebegin + 5,
							mesh.fnodes.data() + _fdist[neighbour[ntria]], mesh.fnodes.data() + _fdist[neighbour[ntria] + 1]);
				}
				continue;
			}

			if (squares == 1 && triangles == 4) {
				size_t ebegin = tnodes.size();
				ttype.push_back((int)Element::CODE::PYRAMID5);
				tsize.push_back(5);
				tnodes.insert(tnodes.end(), 5, -1);
				size_t osquare = obegin, nsquare = nbegin, otria = obegin, ntria = nbegin;
				findElementWithSize(owner, osquare, oindex, 4);
				findElementWithSize(neighbour, nsquare, nindex, 4);
				findElementWithSize(owner, otria, oindex, 3);
				findElementWithSize(neighbour, ntria, nindex, 3);

				if (osquare < oindex) {
					tnodes[ebegin + 0] = mesh.fnodes[_fdist[owner[osquare]] + 0];
					tnodes[ebegin + 1] = mesh.fnodes[_fdist[owner[osquare]] + 1];
					tnodes[ebegin + 2] = mesh.fnodes[_fdist[owner[osquare]] + 2];
					tnodes[ebegin + 3] = mesh.fnodes[_fdist[owner[osquare]] + 3];
				} else {
					tnodes[ebegin + 0] = mesh.fnodes[_fdist[neighbour[nsquare]] + 0];
					tnodes[ebegin + 1] = mesh.fnodes[_fdist[neighbour[nsquare]] + 3];
					tnodes[ebegin + 2] = mesh.fnodes[_fdist[neighbour[nsquare]] + 2];
					tnodes[ebegin + 3] = mesh.fnodes[_fdist[neighbour[nsquare]] + 1];
				}

				if (otria < oindex) {
					tnodes[ebegin + 4] = getUnknown(
							tnodes.data() + ebegin, tnodes.data() + ebegin + 4,
							mesh.fnodes.data() + _fdist[owner[otria]], mesh.fnodes.data() + _fdist[owner[otria] + 1]);
				} else {
					tnodes[ebegin + 4] = getUnknown(
							tnodes.data() + ebegin, tnodes.data() + ebegin + 4,
							mesh.fnodes.data() + _fdist[neighbour[ntria]], mesh.fnodes.data() + _fdist[neighbour[ntria] + 1]);
				}
				continue;
			}
			ESINFO(ERROR) << "OpenFOAM parser: an unknown element type with " << triangles << " triangles and " << squares << " squares [ID=" << e << "].";
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
	mesh.body.resize(tdistribution.back(), 0);
	mesh.material.resize(tdistribution.back(), 0);
	mesh.eIDs.resize(tdistribution.back());
	std::iota(mesh.eIDs.begin(), mesh.eIDs.end(), eoffset);
}

void OpenFOAMLoader::buildFaces(PlainOpenFOAMData &mesh)
{
	for (size_t i = 0; i < _boundaryData.size(); i++) {
		auto &indices = mesh.eregions[_boundaryData[i].name];
		for (size_t f = 0; f < mesh.fID.size(); f++) {
			if (_boundaryData[i].startFace <= mesh.fID[f] && mesh.fID[f] < _boundaryData[i].startFace + _boundaryData[i].nFaces) {
				indices.push_back(_edist.back() + mesh.fID[f]);
				mesh.eIDs.push_back(_edist.back() + mesh.fID[f]);
				mesh.esize.push_back(mesh.fsize[f]);
				if (mesh.fsize[f] == 3) {
					mesh.etype.push_back((int)Element::CODE::TRIANGLE3);
				}
				if (mesh.fsize[f] == 4) {
					mesh.etype.push_back((int)Element::CODE::SQUARE4);
				}
				mesh.enodes.insert(mesh.enodes.end(), mesh.fnodes.begin() + _fdist[f], mesh.fnodes.begin() + _fdist[f + 1]);
			}
		}
	}

	mesh.body.resize(mesh.esize.size(), 0);
	mesh.material.resize(mesh.esize.size(), 0);
}
