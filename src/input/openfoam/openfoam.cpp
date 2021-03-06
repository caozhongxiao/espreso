
#include "openfoam.h"

#include "parser/points.h"
#include "parser/faces.h"
#include "parser/boundary.h"
#include "parser/zones.h"
#include "parser/sets.h"

#include "input/plaindata.h"
#include "basis/containers/tarray.h"
#include "basis/utilities/parser.h"
#include "basis/utilities/utils.h"
#include "esinfo/envinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.hpp"
#include "config/ecf/input/input.h"

#include "mesh/elements/element.h"
#include "input/randominput.h"
#include "input/sequentialinput.h"

#include <numeric>

using namespace espreso;

void OpenFOAMLoader::load(const InputConfiguration &configuration, Mesh &mesh)
{
	OpenFOAMLoader(configuration, mesh);
}

OpenFOAMLoader::OpenFOAMLoader(const InputConfiguration &configuration, Mesh &mesh)
: _configuration(configuration), _loaders(_configuration, *MPITools::procs),
	_points(80), _faces(80), _owner(80), _neighbour(80), _boundary(80),
	_pointZones(80), _faceZones(80), _cellZones(80)
{
	eslog::start("MESIO: LOAD OpenFOAM", "MESIO");
	eslog::param("database", _configuration.path.c_str());
	eslog::ln();

	readData();
	eslog::checkpointln("MESIO: DATA READ");

	PlainOpenFOAMData meshData;
	parseData(meshData);
	eslog::checkpointln("MESIO: DATA PARSED");

	buildFaces(meshData);
	eslog::checkpointln("MESIO: FACE BUILT");

	collectFaces(meshData);
	eslog::checkpointln("MESIO: FACE COLLECTED");

	buildElements(meshData);
	eslog::checkpointln("MESIO: ELEMENTS BUILT");

	if (info::mpi::size > 1) {
		RandomInput::buildMesh(meshData, mesh);
	} else {
		SequentialInput::buildMesh(meshData, mesh);
	}

	eslog::endln("MESIO: MESH BUILT");
}

void OpenFOAMLoader::distributedReader(const std::string &file, ParallelFile &pfile, bool isMandatory)
{
	if (_loaders.within.rank == 0) {
		MPI_File MPIFile;
		if (!MPILoader::open(_loaders.across, MPIFile, _configuration.path + "/constant/polyMesh/" + file)) {
			if (isMandatory) {
				eslog::error("MPI cannot load file '%s/constant/polyMesh/%s'\n", _configuration.path.c_str(), file.c_str());
			}
		} else {
			MPILoader::read(_loaders.across, MPIFile, pfile);
		}
	}

	MPILoader::scatter(_loaders.within, pfile);
	MPILoader::align(*MPITools::procs, pfile, 0);
}

void OpenFOAMLoader::readData()
{
	ProcessesReduction onlyroot;
	onlyroot.granularity = ProcessesReduction::Granularity::PROCESSES;
	onlyroot.pattern = ProcessesReduction::Pattern::PREFIX;
	onlyroot.reduction_ratio = 1;

	MPISubset singleton(onlyroot, *MPITools::procs);

	auto singletonReader = [&] (const std::string &file, ParallelFile &pfile) {
		if (singleton.within.rank == 0) {
			MPI_File MPIFile;
			if (!MPILoader::open(singleton.across, MPIFile, _configuration.path + "/constant/polyMesh/" + file)) {
				eslog::error("MPI cannot load file '%s/constant/polyMesh/%s'\n", _configuration.path.c_str(), file.c_str());
			}
			MPILoader::read(singleton.across, MPIFile, pfile);
		}
		MPILoader::bcast(singleton.within, pfile);
	};

	distributedReader("points", _points, true);
	distributedReader("faces", _faces, true);
	distributedReader("neighbour", _neighbour, true);
	distributedReader("owner", _owner, true);
	singletonReader("boundary", _boundary);

	distributedReader("pointZones", _pointZones, false);
	distributedReader("faceZones", _faceZones, false);
	distributedReader("cellZones", _cellZones, false);

	if (info::mpi::rank == 0) {
		OpenFOAMSets::inspect(_configuration.path + "/constant/polyMesh/sets/*", _sets);
		size_t size = _sets.size();
		MPI_Bcast(&size, sizeof(size_t), MPI_BYTE, 0, info::mpi::comm);
		MPI_Bcast(_sets.data(), _sets.size() * sizeof(OpenFOAMSet), MPI_BYTE, 0, info::mpi::comm);
	} else {
		size_t size;
		MPI_Bcast(&size, sizeof(size_t), MPI_BYTE, 0, info::mpi::comm);
		_sets.resize(size);
		MPI_Bcast(_sets.data(), _sets.size() * sizeof(OpenFOAMSet), MPI_BYTE, 0, info::mpi::comm);
	}
}

void OpenFOAMLoader::parseData(PlainOpenFOAMData &mesh)
{
	if (!OpenFOAMPoints(_points.begin, _points.end).readData(mesh.nIDs, mesh.coordinates, _configuration.scale_factor)) {
		eslog::error("OpenFOAM loader: cannot parse points.\n");
	}
	if (!OpenFOAMFaces(_faces.begin, _faces.end).readFaces(mesh)) {
		eslog::error("OpenFOAM loader: cannot parse faces.\n");
	}
	if (!OpenFOAMFaces(_owner.begin, _owner.end).readParents(mesh.owner)) {
		eslog::error("OpenFOAM loader: cannot parse owners.\n");
	}
	if (!OpenFOAMFaces(_neighbour.begin, _neighbour.end).readParents(mesh.neighbour)) {
		eslog::error("OpenFOAM loader: cannot parse neighbours.\n");
	}

	if (!OpenFOAMBoundary(_boundary.begin, _boundary.end).readData(mesh)) {
		eslog::error("OpenFOAM loader: cannot parse boundary.\n");
	}

	if (_pointZones.offsets.back() != 0 && !OpenFOAMZones(_pointZones).readPoints(mesh)) {
		eslog::error("OpenFOAM loader: cannot parse pointZones.\n");
	}

	if (_faceZones.offsets.back() != 0 && !OpenFOAMZones(_faceZones).readFaces(mesh)) {
		eslog::error("OpenFOAM loader: cannot parse faceZones.\n");
	}

	if (_cellZones.offsets.back() != 0 && !OpenFOAMZones(_cellZones).readCells(mesh)) {
		eslog::error("OpenFOAM loader: cannot parse cellZones.\n");
	}
}

void OpenFOAMLoader::buildFaces(PlainOpenFOAMData &mesh)
{
	// 1. Find MAX element ID in order to be able correctly set face IDs in continuous interval

	esint maxID = 0;
	if (mesh.owner.size()) {
		maxID = *std::max_element(mesh.owner.begin(), mesh.owner.end());
	}
	if (mesh.neighbour.size()) {
		maxID = std::max(*std::max_element(mesh.neighbour.begin(), mesh.neighbour.end()), maxID);
	}

	MPI_Allreduce(&maxID, &mesh.nelements, sizeof(esint), MPI_BYTE, MPITools::esintOperations().max, info::mpi::comm);
	mesh.nelements += 1;

	_edist = tarray<esint>::distribute(info::mpi::size, mesh.nelements);

	std::vector<esint> fDistribution = Communication::getDistribution<esint>(mesh.fsize.size());

	// 2. Add sets that are not in any zone or boundary

	for (size_t i = 0; i < _sets.size(); i++) {
		std::string name = std::string(_sets[i].name);
		switch (_sets[i].type) {
		case OpenFOAMSet::SetType::CELL_SET:
			for (auto ereg = mesh.eregions.begin(); ereg != mesh.eregions.end(); ++ereg) {
				if (StringCompare::caseInsensitivePreffix(ereg->first, mesh.elementprefix + name)) {
					_sets.erase(_sets.begin() + i--);
					break;
				}
			}
			break;
		case OpenFOAMSet::SetType::FACE_SET:
			for (auto ereg = mesh.eregions.begin(); ereg != mesh.eregions.end(); ++ereg) {
				if (StringCompare::caseInsensitivePreffix(ereg->first, mesh.boundaryprefix + name)) {
					_sets.erase(_sets.begin() + i--);
					break;
				}
			}
			break;
		case OpenFOAMSet::SetType::POINT_SET:
			for (auto nreg = mesh.nregions.begin(); nreg != mesh.nregions.end(); ++nreg) {
				if (StringCompare::caseInsensitivePreffix(nreg->first, name)) {
					_sets.erase(_sets.begin() + i--);
					break;
				}
			}
			break;
		}
	}

	for (size_t i = 0; i < _sets.size(); i++) {
		std::string name = std::string(_sets[i].name);
		switch (_sets[i].type) {
		case OpenFOAMSet::SetType::CELL_SET: {
			ParallelFile file(80);
			distributedReader("/sets/" + name, file, true);
			if (!OpenFOAMSets(file.begin, file.end).readData(_sets[i], mesh.eregions[mesh.elementprefix + name])) {
				eslog::error("OpenFOAM loader: cannot parse set '%s'.\n", name.c_str());
			}
		} break;
		case OpenFOAMSet::SetType::FACE_SET: {
			ParallelFile file(80);
			distributedReader("/sets/" + name, file, true);
			if (!OpenFOAMSets(file.begin, file.end).readData(_sets[i], mesh.eregions[mesh.boundaryprefix + name])) {
				eslog::error("OpenFOAM loader: cannot parse set '%s'.\n", name.c_str());
			}
		} break;
		case OpenFOAMSet::SetType::POINT_SET: {
			ParallelFile file(80);
			distributedReader("/sets/" + name, file, true);
			if (!OpenFOAMSets(file.begin, file.end).readData(_sets[i], mesh.nregions[name])) {
				eslog::error("OpenFOAM loader: cannot parse set '%s'.\n", name.c_str());
			}
		} break;
		}
	}

	// 3. Exchange region data to processes that hold given faces

	std::vector<esint> sBuffer, rBuffer;

	std::vector<size_t> rpointer(mesh.eregions.size());
	for (int r = 0; r < info::mpi::size; r++) {
		size_t prevsize = sBuffer.size();
		sBuffer.push_back(0); // total size
		sBuffer.push_back(r); // target

		size_t rindex = 0;
		for (auto it = mesh.eregions.begin(); it != mesh.eregions.end(); ++it, ++rindex) {
			if (StringCompare::caseInsensitivePreffix(mesh.boundaryprefix, it->first)) {
				size_t prevrsize = sBuffer.size();
				sBuffer.push_back(0); // region size

				for ( ; rpointer[rindex] < it->second.size() && it->second[rpointer[rindex]] < fDistribution[r + 1]; ++rpointer[rindex]) {
					sBuffer.push_back(it->second[rpointer[rindex]]);
				}
				sBuffer[prevrsize] = sBuffer.size() - prevrsize - 1;
			}
		}

		sBuffer[prevsize] = sBuffer.size() - prevsize;
	}

	if (!Communication::allToAllWithDataSizeAndTarget(sBuffer, rBuffer)) {
		eslog::error("ESPRESO internal error: distribute faces indices in regions.\n");
	}

	for (auto it = mesh.eregions.begin(); it != mesh.eregions.end(); ++it) {
		if (StringCompare::caseInsensitivePreffix(mesh.boundaryprefix, it->first)) {
			it->second.clear();
		}
	}

	std::vector<esint> usedfaces;
	size_t offset = 0;
	for (int r = 0; r < info::mpi::size; r++) {
		++offset; // skip total size
		++offset; // skip target

		size_t rindex = 0;
		for (auto it = mesh.eregions.begin(); it != mesh.eregions.end(); ++it, ++rindex) {
			if (StringCompare::caseInsensitivePreffix(mesh.boundaryprefix, it->first)) {
				size_t rsize = rBuffer[offset++];

				for (size_t i = 0; i < rsize; ++i) {
					usedfaces.push_back(rBuffer[offset]);
					it->second.push_back(rBuffer[offset++]); // index faces after elements
				}
			}
			std::sort(it->second.begin(), it->second.end());
		}
	}

	// 4. Add used faces into elements

	utils::sortAndRemoveDuplicity(usedfaces);

	size_t foffset = usedfaces.size();
	Communication::exscan(foffset);
	foffset += mesh.nelements;

	_fdist.reserve(mesh.fsize.size() + 1);
	_fdist.push_back(0);
	for (size_t f = 0; f < mesh.fsize.size(); f++) {
		_fdist.push_back(_fdist.back() + mesh.fsize[f]);
	}

	for (size_t i = 0; i < usedfaces.size(); i++) {
		esint findex = usedfaces[i] - fDistribution[info::mpi::rank];

		mesh.esize.push_back(mesh.fsize[findex]);
		mesh.enodes.insert(mesh.enodes.end(), mesh.fnodes.begin() + _fdist[findex], mesh.fnodes.begin() + _fdist[findex + 1]);
		if (mesh.fsize[findex] == 3) {
			mesh.etype.push_back((int)Element::CODE::TRIANGLE3);
		}
		if (mesh.fsize[findex] == 4) {
			mesh.etype.push_back((int)Element::CODE::SQUARE4);
		}
	}

	mesh.eIDs.resize(mesh.esize.size(), 0);
	std::iota(mesh.eIDs.begin(), mesh.eIDs.end(), foffset);
	mesh.body.resize(mesh.esize.size(), 0);
	mesh.material.resize(mesh.esize.size(), 0);

	for (auto it = mesh.eregions.begin(); it != mesh.eregions.end(); ++it) {
		if (StringCompare::caseInsensitivePreffix(mesh.boundaryprefix, it->first)) {
			auto id = usedfaces.begin();
			for (size_t i = 0; i < it->second.size(); i++) {
				while (*id < it->second[i]) { ++id; }
				it->second[i] = foffset + id - usedfaces.begin();
			}
		}
	}
}

void OpenFOAMLoader::collectFaces(PlainOpenFOAMData &mesh)
{
	std::vector<size_t> ownersDist = Communication::getDistribution(mesh.owner.size());
	std::vector<size_t> neighborsDist = Communication::getDistribution(mesh.neighbour.size());
	std::vector<size_t> target = Communication::getDistribution(mesh.fsize.size());

	if (!Communication::balance(mesh.owner, ownersDist, target)) {
		eslog::error("ESPRESO internal error: balance faces owners.\n");
	}

	for (size_t i = 0; i < target.size(); i++) {
		if (target[i] > neighborsDist.back()) {
			target[i] = neighborsDist.back();
		}
	}

	if (!Communication::balance(mesh.neighbour, neighborsDist, target)) {
		eslog::error("ESPRESO internal error: balance faces neighbors.\n");
	}

	size_t firstID = Communication::getDistribution(mesh.fsize.size())[info::mpi::rank];

	auto sortIDs = [&] (std::vector<esint> &permutation, const std::vector<esint> &data) {
		permutation.resize(data.size());
		std::iota(permutation.begin(), permutation.end(), 0);

		std::sort(permutation.begin(), permutation.end(), [&] (esint i, esint j) {
			if (data[i] == data[j]) {
				return i < j;
			}
			return data[i] < data[j];
		});
	};

	std::vector<esint> oPermutation, nPermutation;
	sortIDs(oPermutation, mesh.owner);
	sortIDs(nPermutation, mesh.neighbour);

	std::vector<esint> sBuffer, rBuffer;
	// ID, size, owner / -1 * neighbor, nodes
	sBuffer.reserve(4 * info::mpi::size + 3 * (mesh.owner.size() + mesh.neighbour.size()) + mesh.fnodes.size());

	size_t prevsize;
	auto obegin = oPermutation.begin();
	auto nbegin = nPermutation.begin();
	for (int r = 0; r < info::mpi::size; r++) {
		prevsize = sBuffer.size();
		sBuffer.push_back(0); // total size
		sBuffer.push_back(r); // target
		sBuffer.push_back(0); // number of faces

		auto o = obegin;
		for ( ; o != oPermutation.end() && mesh.owner[*o] < _edist[r + 1]; ++o) {
			sBuffer.push_back(firstID + *o);
			sBuffer.push_back(mesh.owner[*o]);
			sBuffer.push_back(mesh.fsize[*o]);
			sBuffer.insert(sBuffer.end(), mesh.fnodes.begin() + _fdist[*o], mesh.fnodes.begin() + _fdist[*o + 1]);
		}
		auto n = nbegin;
		for ( ; n != nPermutation.end() && mesh.neighbour[*n] < _edist[r + 1]; ++n) {
			sBuffer.push_back(firstID + *n);
			sBuffer.push_back(-mesh.neighbour[*n] - 1);
			sBuffer.push_back(mesh.fsize[*n]);
			sBuffer.insert(sBuffer.end(), mesh.fnodes.begin() + _fdist[*n], mesh.fnodes.begin() + _fdist[*n + 1]);
		}
		sBuffer[prevsize + 2] = (o - obegin) + (n - nbegin);
		obegin = o;
		nbegin = n;

		sBuffer[prevsize] = sBuffer.size() - prevsize;
	}

	if (!Communication::allToAllWithDataSizeAndTarget(sBuffer, rBuffer)) {
		eslog::error("ESPRESO internal error: distribute permuted elements.\n");
	}

	mesh.fIDs.clear();
	mesh.fsize.clear();
	mesh.fnodes.clear();
	mesh.owner.clear();
	mesh.neighbour.clear();

	std::vector<esint> fIDs, fsize, fnodes, owners;

	size_t offset = 0;
	for (int r = 0; r < info::mpi::size; r++) {
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

	std::vector<esint> fpermutation(fIDs.size());
	std::iota(fpermutation.begin(), fpermutation.end(), 0);
	std::sort(fpermutation.begin(), fpermutation.end(), [&] (esint i, esint j) {
		return fIDs[i] < fIDs[j];
	});

	_fdist.resize(1);
	_fdist.reserve(fsize.size() + 1);
	for (size_t e = 0; e < fsize.size(); e++) {
		_fdist.push_back(_fdist.back() + fsize[e]);
	}

	if (fpermutation.size()) {
		size_t rest = 0;
		for (auto i = fpermutation.begin(); i != fpermutation.end() - 1; ++i) {
			if (fIDs[*i] == fIDs[*(i + 1)]) {
				mesh.fsize.push_back(fsize[*i]);
				if (owners[*i] >= 0) {
					mesh.owner.push_back(owners[*i]);
					mesh.neighbour.push_back(-owners[*(i + 1)] - 1);
					mesh.fnodes.insert(mesh.fnodes.end(), fnodes.begin() + _fdist[*i], fnodes.begin() + _fdist[*i + 1]);
				} else {
					mesh.owner.push_back(owners[*(i + 1)]);
					mesh.neighbour.push_back(-owners[*i] - 1);
					mesh.fnodes.insert(mesh.fnodes.end(), fnodes.begin() + _fdist[*(i + 1)], fnodes.begin() + _fdist[*(i + 1) + 1]);
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
			if (owners[*i] >= 0) {
				mesh.owner.push_back(owners[*i]);
				mesh.fnodes.insert(mesh.fnodes.end(), fnodes.begin() + _fdist[*i], fnodes.begin() + _fdist[*i + 1]);
			} else {
				mesh.owner.push_back(-owners[*i] - 1);
				mesh.fnodes.insert(mesh.fnodes.end(), fnodes.rbegin() + _fdist.back() - _fdist[*i + 1], fnodes.rbegin() + _fdist.back() - _fdist[*i]);
			}
		}
	}
}

void OpenFOAMLoader::buildElements(PlainOpenFOAMData &mesh)
{
	size_t threads = info::env::OMP_NUM_THREADS;

	auto sortIDs = [&] (std::vector<esint> &permutation, const std::vector<esint> &data) {
		permutation.resize(data.size());
		std::iota(permutation.begin(), permutation.end(), 0);

		std::sort(permutation.begin(), permutation.end(), [&] (esint i, esint j) {
			if (data[i] == data[j]) {
				return i < j;
			}
			return data[i] < data[j];
		});
	};

	std::vector<esint> owner, neighbour;
	sortIDs(owner, mesh.owner);
	sortIDs(neighbour, mesh.neighbour);

	std::vector<size_t> tdistribution = tarray<size_t>::distribute(threads, _edist[info::mpi::rank + 1] - _edist[info::mpi::rank]);

	std::vector<std::vector<esint> > esize(threads), enodes(threads);
	std::vector<std::vector<int> > etype(threads);

	_fdist.clear();
	_fdist.reserve(tdistribution.back() + 1);
	_fdist.push_back(0);
	for (size_t f = 0; f < mesh.fsize.size(); f++) {
		_fdist.push_back(_fdist.back() + mesh.fsize[f]);
	}

	auto getThreadBegin = [] (const std::vector<esint> &data, const std::vector<esint> &perm, esint eindex) {
		return std::lower_bound(perm.begin(), perm.end(), eindex, [&] (esint i, esint eindex) {
			return data[i] < eindex;
		}) - perm.begin();
	};

	auto addFaces = [&] (const std::vector<esint> &data, const std::vector<esint> &perm, size_t &triangles, size_t &squares, size_t &index, esint element) {
		while (index < perm.size() && data[perm[index]] == element) {
			switch (mesh.fsize[perm[index++]]) {
			case 3: ++triangles; break;
			case 4:   ++squares; break;
			}
		}
	};

	auto getFace = [&] (const std::vector<esint> &data, const std::vector<esint> &perm, size_t index, esint element, esint n1, esint n2) {
		while (index < perm.size() && data[perm[index]] == element) {
			for (esint f = 0; f < mesh.fsize[perm[index]]; f++) {
				if (n1 == mesh.fnodes[_fdist[perm[index]] + f] && n2 == mesh.fnodes[_fdist[perm[index]] + (f + 1) % mesh.fsize[perm[index]]]) {
					return std::pair<size_t, esint>(index, f);
				}
			}
			++index;
		}
		return std::pair<size_t, esint>(index, -1);
	};

	auto getUnknown = [&] (esint *kbegin, esint *kend, esint *ubegin, esint *uend) {
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
		return (esint)-1;
	};

	auto findElementWithSize = [&] (const std::vector<esint> &perm, size_t &index, size_t &max, int size) {
		while (index < max) { // there is at least one owner
			if (mesh.fsize[perm[index]] == size) {
				break;
			} else {
				++index;
			}
		}
	};

	size_t eoffset = _edist[info::mpi::rank];

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<esint> tsize, tnodes;
		std::vector<int> ttype;

		size_t oindex = getThreadBegin(mesh.owner, owner, tdistribution[t] + eoffset);
		size_t nindex = getThreadBegin(mesh.neighbour, neighbour, tdistribution[t] + eoffset);

		for (size_t e = tdistribution[t] + eoffset; e < tdistribution[t + 1] + eoffset; e++) {
			size_t obegin = oindex, nbegin = nindex;
			std::pair<size_t, esint> index;
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
			eslog::error("OpenFOAM parser: an unknown element type with '%ld' triangles and '%ld' squares [ID='%ld'].\n", triangles, squares, e);
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
	mesh.body.resize(mesh.esize.size(), 0);
	mesh.material.resize(mesh.esize.size(), 0);
	size_t fsize = mesh.eIDs.size();
	mesh.eIDs.resize(mesh.esize.size());
	std::iota(mesh.eIDs.begin() + fsize, mesh.eIDs.end(), eoffset);
}


