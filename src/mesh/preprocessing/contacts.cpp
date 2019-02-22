
#include "meshpreprocessing.h"

#include "mesh/mesh.h"
#include "mesh/elements/element.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/surfacestore.h"
#include "mesh/store/contactstore.h"

#include "basis/containers/serializededata.h"
#include "basis/matrices/denseMatrix.h"
#include "basis/utilities/communication.h"
#include "basis/utilities/utils.h"

#include "esinfo/mpiinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.h"

#include <algorithm>
#include <numeric>

#include "basis/utilities/print.h"

using namespace espreso;


void MeshPreprocessing::computeBodiesSurface()
{
	if (_mesh->elements->neighbors == NULL) {
		this->computeElementsNeighbors();
	}

	eslog::startln("MESH: COMPUTE SURFACE");

	size_t threads = info::env::OMP_NUM_THREADS;

	std::vector<std::vector<esint> > faces(threads), facesDistribution(threads), ecounters(threads, std::vector<esint>((int)Element::CODE::SIZE));
	std::vector<std::vector<Element*> > fpointers(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto nodes = _mesh->elements->procNodes->cbegin(t);
		auto neighs = _mesh->elements->neighbors->cbegin(t);
		const auto &epointers = _mesh->elements->epointers->datatarray();

		std::vector<esint> fdist, fdata, ecounter((int)Element::CODE::SIZE);
		std::vector<Element*> fpointer;
		if (t == 0) {
			fdist.push_back(0);
		}

		for (size_t e = _mesh->elements->distribution[t]; e < _mesh->elements->distribution[t + 1]; ++e, ++neighs, ++nodes) {
			for (size_t n = 0; n < neighs->size(); ++n) {
				if (neighs->at(n) == -1) {
					auto face = epointers[e]->faces->begin() + n;
					for (auto f = face->begin(); f != face->end(); ++f) {
						fdata.push_back(nodes->at(*f));
					}
					fdist.push_back(fdata.size());
					fpointer.push_back(epointers[e]->facepointers->datatarray()[n]);
					++ecounter[(int)fpointer.back()->code];
				}
			}
		}

		facesDistribution[t].swap(fdist);
		faces[t].swap(fdata);
		fpointers[t].swap(fpointer);
		ecounters[t].swap(ecounter);
	}

	if (_mesh->surface == NULL) {
		_mesh->surface = new SurfaceStore();
	}

	for (size_t t = 1; t < threads; t++) {
		for (size_t e = 0; e < ecounters[0].size(); e++) {
			ecounters[0][e] += ecounters[t][e];
		}
	}

	serializededata<esint, Element*>::balance(1, fpointers);
	_mesh->surface->epointers = new serializededata<esint, Element*>(1, fpointers);
	_mesh->surface->ecounters = ecounters[0];

	_mesh->surface->edistribution = _mesh->surface->epointers->datatarray().distribution();

	if (_mesh->surface->ecounters[(int)Element::CODE::TRIANGLE3] == (esint)_mesh->surface->edistribution.back()) {
		serializededata<esint, esint>::balance(3, faces, &_mesh->surface->edistribution);
		_mesh->surface->elements = new serializededata<esint, esint>(3, faces);
		_mesh->surface->triangles = _mesh->surface->elements;
		_mesh->surface->tdistribution = _mesh->surface->edistribution;
	} else {
		utils::threadDistributionToFullDistribution(facesDistribution);
		serializededata<esint, esint>::balance(facesDistribution, faces, &_mesh->surface->edistribution);
		_mesh->surface->elements = new serializededata<esint, esint>(facesDistribution, faces);
	}

	eslog::endln("MESH: SURFACE COMPUTED");
}

void MeshPreprocessing::computeSurfaceLocations()
{
	eslog::startln("MESH: COMPUTE SURFACE LOCATIONS");

	size_t threads = info::env::OMP_NUM_THREADS;

	std::vector<std::vector<Point> > points(threads);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		points[t].resize(_mesh->surface->elements->cbegin(t + 1)->begin() - _mesh->surface->elements->cbegin(t)->begin());
		size_t i = 0;
		for (auto e = _mesh->surface->elements->cbegin(t); e != _mesh->surface->elements->cend(t); ++e) {
			for (auto n = e->begin(); n != e->end(); ++n, ++i) {
				points[t][i] = _mesh->nodes->coordinates->datatarray()[*n];
			}
		}
	}

	if (_mesh->surface->elements->boundarytarray().size()) {
		_mesh->contacts->elements = new serializededata<esint, Point>(tarray<esint>(_mesh->surface->elements->boundarytarray()), points);
	} else {
		_mesh->contacts->elements = new serializededata<esint, Point>(3, points);
	}

	size_t precision = 0;
	double epsilon = 1e-6;
	while (true) {
		double value = _mesh->contacts->eps * pow(10, precision);
		if (std::round(value) <= value + epsilon && value - epsilon <= std::round(value)) {
			break;
		} else {
			++precision;
		}
	}

	_mesh->contacts->boundingBox[0] = Point();
	_mesh->contacts->boundingBox[1] = Point();

	std::vector<Point> bbox(2 * _mesh->nodes->externalIntervals.size());

	// 1. compute bounding box
	//////////////////////////

	#pragma omp parallel for
	for (size_t t = 0; t < _mesh->nodes->externalIntervals.size(); t++) {
		esint eint = _mesh->nodes->externalIntervals[t];
		Point min = _mesh->nodes->coordinates->datatarray().front();
		Point max = _mesh->nodes->coordinates->datatarray().front();
		for (esint n = _mesh->nodes->pintervals[eint].begin; n < _mesh->nodes->pintervals[eint].end; ++n) {
			min.x = std::min(min.x, _mesh->nodes->coordinates->datatarray()[n].x);
			min.y = std::min(min.y, _mesh->nodes->coordinates->datatarray()[n].y);
			min.z = std::min(min.z, _mesh->nodes->coordinates->datatarray()[n].z);
			max.x = std::max(max.x, _mesh->nodes->coordinates->datatarray()[n].x);
			max.y = std::max(max.y, _mesh->nodes->coordinates->datatarray()[n].y);
			max.z = std::max(max.z, _mesh->nodes->coordinates->datatarray()[n].z);
		}

		bbox[2 * t + 0] = min;
		bbox[2 * t + 1] = max;
	}

	for (size_t t = 1; t < _mesh->nodes->externalIntervals.size(); t++) {
		bbox[0].x = std::min(bbox[0].x, bbox[2 * t + 0].x);
		bbox[0].y = std::min(bbox[0].y, bbox[2 * t + 0].y);
		bbox[0].z = std::min(bbox[0].z, bbox[2 * t + 0].z);
		bbox[1].x = std::max(bbox[1].x, bbox[2 * t + 1].x);
		bbox[1].y = std::max(bbox[1].y, bbox[2 * t + 1].y);
		bbox[1].z = std::max(bbox[1].z, bbox[2 * t + 1].z);
	}

	if (_mesh->nodes->externalIntervals.size()) {
		_mesh->contacts->boundingBox[0] = bbox[0];
		_mesh->contacts->boundingBox[1] = bbox[1];
	}

	auto rounddown = [&] (double &value) {
		int rounder = _mesh->contacts->eps * std::pow(10, precision);
		int result = std::ceil(value * pow(10, precision));
		result = result - (rounder - result % rounder);
		value = result / (double)std::pow(10, precision);
	};
	auto roundup = [&] (double &value) {
		int rounder = _mesh->contacts->eps * std::pow(10, precision);
		int result = std::floor(value * pow(10, precision));
		result = result + (rounder - result % rounder);
		value = result / (double)std::pow(10, precision);
	};

	rounddown(_mesh->contacts->boundingBox[0].x);
	rounddown(_mesh->contacts->boundingBox[0].y);
	rounddown(_mesh->contacts->boundingBox[0].z);
	roundup(_mesh->contacts->boundingBox[1].x);
	roundup(_mesh->contacts->boundingBox[1].y);
	roundup(_mesh->contacts->boundingBox[1].z);

	size_t xmaxsize = std::ceil((_mesh->contacts->boundingBox[1].x - _mesh->contacts->boundingBox[0].x) / _mesh->contacts->eps);
	size_t ymaxsize = std::ceil((_mesh->contacts->boundingBox[1].y - _mesh->contacts->boundingBox[0].y) / _mesh->contacts->eps);
	size_t zmaxsize = std::ceil((_mesh->contacts->boundingBox[1].z - _mesh->contacts->boundingBox[0].z) / _mesh->contacts->eps);

	double avgsize = (xmaxsize + ymaxsize + zmaxsize) / 3.;
	_mesh->contacts->groupsize = std::max(avgsize / std::pow(_mesh->contacts->surface->elements->structures(), 1. / 3) / 1.5, 1.);

	MPI_Allreduce(&_mesh->contacts->boundingBox[0], &_mesh->contacts->globalBox[0], 3, MPI_DOUBLE, MPI_MIN, info::mpi::comm);
	MPI_Allreduce(&_mesh->contacts->boundingBox[1], &_mesh->contacts->globalBox[1], 3, MPI_DOUBLE, MPI_MAX, info::mpi::comm);

	std::vector<std::vector<std::pair<esint, esint> > > grids(threads);

	double boxsize = _mesh->contacts->eps * _mesh->contacts->groupsize;
	size_t xsize = std::ceil((_mesh->contacts->boundingBox[1].x - _mesh->contacts->boundingBox[0].x) / boxsize);
	size_t ysize = std::ceil((_mesh->contacts->boundingBox[1].y - _mesh->contacts->boundingBox[0].y) / boxsize);

	// 2. map into grid
	///////////////////
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<std::pair<esint, esint> > grid;

		int xmin, xmax, ymin, ymax, zmin, zmax;
		size_t gsize = 0;

		auto insert = [&] (Point &min, Point &max, esint eindex) {
			xmin = std::floor((min.x - _mesh->contacts->boundingBox[0].x - epsilon) / boxsize);
			xmax = std::ceil((max.x - _mesh->contacts->boundingBox[0].x + epsilon) / boxsize);
			ymin = std::floor((min.y - _mesh->contacts->boundingBox[0].y - epsilon) / boxsize);
			ymax = std::ceil((max.y - _mesh->contacts->boundingBox[0].y + epsilon) / boxsize);
			zmin = std::floor((min.z - _mesh->contacts->boundingBox[0].z - epsilon) / boxsize);
			zmax = std::ceil((max.z - _mesh->contacts->boundingBox[0].z + epsilon) / boxsize);

			grid.resize(gsize + (xmax - xmin) * (ymax - ymin) * (zmax - zmin), {0, eindex});
			for (int z = zmin; z < zmax; ++z) {
				for (int y = ymin; y < ymax; ++y) {
					for (int x = xmin; x < xmax; ++x, ++gsize) {
						grid[gsize].first = xsize * ysize * z + xsize * y + x;
					}
				}
			}
		};

		Point min, max;
		esint eindex = _mesh->contacts->surface->edistribution[t];
		for (auto e = _mesh->contacts->elements->cbegin(t); e != _mesh->contacts->elements->cend(t); ++e, ++eindex) {
			min = max = e->front();
			for (auto n = e->begin() + 1; n != e->end(); ++n) {
				min.x = std::min(min.x, n->x);
				min.y = std::min(min.y, n->y);
				min.z = std::min(min.z, n->z);
				max.x = std::max(max.x, n->x);
				max.y = std::max(max.y, n->y);
				max.z = std::max(max.z, n->z);
			}
			insert(min, max, eindex);
		}

		grids[t].swap(grid);
	}

	{
		std::vector<size_t> gdistribution = { 0, grids[0].size() };
		for (size_t t = 1; t < threads; t++) {
			grids[0].insert(grids[0].end(), grids[t].begin(), grids[t].end());
			gdistribution.push_back(grids[0].size());
		}
		utils::sortWithInplaceMerge(grids[0], gdistribution);
	}

	// 3. create sparse structure
	///////////////////

	std::vector<size_t> distribution = tarray<size_t>::distribute(threads, grids[0].size());

	for (size_t t = 1; t < threads; t++) {
		while (distribution[t] != distribution[t + 1] && grids[0][distribution[t]].first == grids[0][distribution[t] - 1].first) {
			++distribution[t];
		}
	}

	std::vector<std::vector<esint> > gdist(threads), gdata(threads), gfilled(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<esint> dist, data, filled;
		if (t == 0) {
			dist.push_back(0);
		}
		if (distribution[t] != distribution[t + 1]) {
			filled.push_back(grids[0][distribution[t]].first);
		}

		for (size_t i = distribution[t]; i < distribution[t + 1]; ++i) {
			if (filled.back() != grids[0][i].first) {
				filled.push_back(grids[0][i].first);
				dist.push_back(data.size());
			}
			data.push_back(grids[0][i].second);
		}
		dist.push_back(data.size());

		gdist[t].swap(dist);
		gdata[t].swap(data);
		gfilled[t].swap(filled);
	}

	utils::threadDistributionToFullDistribution(gdist);
	for (size_t t = 1; t < threads; t++) {
		gdist[0].insert(gdist[0].end(), gdist[t].begin(), gdist[t].end());
		gdata[0].insert(gdata[0].end(), gdata[t].begin(), gdata[t].end());
		gfilled[0].insert(gfilled[0].end(), gfilled[t].begin(), gfilled[t].end());
	}

	// never threaded
	gdist.resize(1);
	gdata.resize(1);

	_mesh->contacts->filledCells.swap(gfilled[0]);
	_mesh->contacts->grid = new serializededata<esint, esint>(gdist, gdata);

	eslog::endln("MESH: SURFACE LOCATIONS COMPUTED");
}

void MeshPreprocessing::computeContactNormals()
{
	if (_mesh->contacts->elements == NULL) {
		_mesh->preprocessing->computeSurfaceLocations();
	}

	eslog::startln("MESH: COMPUTE CONTACT NORMALS");

	size_t threads = info::env::OMP_NUM_THREADS;

	if (_mesh->contacts->enormals != NULL) {
		delete _mesh->contacts->enormals;
	}
	_mesh->contacts->enormals = new serializededata<esint, Point>(*_mesh->contacts->elements);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		const auto &epointers = _mesh->surface->epointers->datatarray();
		auto element = _mesh->contacts->elements->begin(t);
		auto normals = _mesh->contacts->enormals->begin(t);
		DenseMatrix tangents(2, 3), coords(8, 3); // maximal size

		for (auto e = _mesh->contacts->surface->edistribution[t]; e < _mesh->contacts->surface->edistribution[t + 1]; ++e, ++element, ++normals) {
			coords.resize(epointers[e]->nodes, 3);
			int cindex = 0;
			for (auto n = element->begin(); n != element->end(); ++n, ++cindex) {
				coords(cindex, 0) = n->x;
				coords(cindex, 1) = n->y;
				coords(cindex, 2) = n->z;
			}
			for (int n = 0; n < epointers[e]->nodes; n++) {
				tangents.multiply((*epointers[e]->ndN)[n], coords);
				Point t(tangents(0, 0), tangents(0, 1), tangents(0, 2));
				Point s(tangents(1, 0), tangents(1, 1), tangents(1, 2));
				normals->at(n) = Point::cross(t, s).normalize();
			}
		}
	}

	eslog::endln("MESH: CONTACT NORMALS COMPUTED");
}

void MeshPreprocessing::searchContactInterfaces()
{
	eslog::startln("MESH: SEARCH CONTACT INTERFACE");

	size_t threads = info::env::OMP_NUM_THREADS;
	double epsilon = 1e-6;

	double boxsize = _mesh->contacts->eps * _mesh->contacts->groupsize;
	int xsize = std::ceil((_mesh->contacts->boundingBox[1].x - _mesh->contacts->boundingBox[0].x) / boxsize);
	int ysize = std::ceil((_mesh->contacts->boundingBox[1].y - _mesh->contacts->boundingBox[0].y) / boxsize);
	int zsize = std::ceil((_mesh->contacts->boundingBox[1].z - _mesh->contacts->boundingBox[0].z) / boxsize);

	auto areIntersected = [&] (Point *block1, Point *block2) {
		return	!
				(block1[1].x + epsilon < block2[0].x || block2[1].x + epsilon < block1[0].x) ||
				(block1[1].y + epsilon < block2[0].y || block2[1].y + epsilon < block1[0].y) ||
				(block1[1].z + epsilon < block2[0].z || block2[1].z + epsilon < block1[0].z);
	};

	auto intersect = [&] (Point *intersection, Point *block1, Point *block2) {
		intersection[0].x = std::max(block1[0].x, block2[0].x);
		intersection[0].y = std::max(block1[0].y, block2[0].y);
		intersection[0].z = std::max(block1[0].z, block2[0].z);

		intersection[1].x = std::min(block1[1].x, block2[1].x);
		intersection[1].y = std::min(block1[1].y, block2[1].y);
		intersection[1].z = std::min(block1[1].z, block2[1].z);
	};

	std::vector<Point> boxes(2 * info::mpi::size);

	MPI_Allgather(_mesh->contacts->boundingBox, 6, MPI_DOUBLE, boxes.data(), 6, MPI_DOUBLE, info::mpi::comm);

	for (int r = 0; r < info::mpi::size; r++) {
		if (r != info::mpi::rank) {
			if (areIntersected(_mesh->contacts->boundingBox, boxes.data() + 2 * r)) {
				_mesh->contacts->gneighbors.push_back(r);
			}
		}
	}

	// EXCHANGE GROUP SIZES
	std::vector<std::vector<size_t> > sGroupSize(_mesh->contacts->gneighbors.size(), { _mesh->contacts->groupsize }), rGroupSize(_mesh->contacts->gneighbors.size(), { 0 });
	if (!Communication::exchangeKnownSize(sGroupSize, rGroupSize, _mesh->contacts->gneighbors)) {
		eslog::error("ESPRESO internal error: exchange group size.\n");
	}

	std::vector<std::vector<esint> > sFilled(_mesh->contacts->gneighbors.size());
	std::vector<std::vector<esint> > sBlock(_mesh->contacts->gneighbors.size()), rBlock(_mesh->contacts->gneighbors.size());
	std::vector<std::vector<esint> > sIDs(_mesh->contacts->gneighbors.size()), rIDs(_mesh->contacts->gneighbors.size());
	std::vector<std::vector<esint> > sDist(_mesh->contacts->gneighbors.size()), rDist(_mesh->contacts->gneighbors.size());
	std::vector<std::vector<Point> > sData(_mesh->contacts->gneighbors.size()), rData(_mesh->contacts->gneighbors.size());

	for (size_t n = 0; n < _mesh->contacts->gneighbors.size(); n++) {
		Point intersection[2];
		intersect(intersection, _mesh->contacts->boundingBox, boxes.data() + 2 * _mesh->contacts->gneighbors[n]);

		int xmin = std::max((int)std::floor((intersection[0].x - _mesh->contacts->boundingBox[0].x - epsilon) / boxsize), 0);
		int ymin = std::max((int)std::floor((intersection[0].y - _mesh->contacts->boundingBox[0].y - epsilon) / boxsize), 0);
		int zmin = std::max((int)std::floor((intersection[0].z - _mesh->contacts->boundingBox[0].z - epsilon) / boxsize), 0);
		int xmax = std::min((int)std::ceil ((intersection[1].x - _mesh->contacts->boundingBox[0].x + epsilon) / boxsize), xsize);
		int ymax = std::min((int)std::ceil ((intersection[1].y - _mesh->contacts->boundingBox[0].y + epsilon) / boxsize), xsize);
		int zmax = std::min((int)std::ceil ((intersection[1].z - _mesh->contacts->boundingBox[0].z + epsilon) / boxsize), xsize);

		auto begin = std::lower_bound(_mesh->contacts->filledCells.begin(), _mesh->contacts->filledCells.end(), xsize * ysize * zmin);
		auto end   = std::lower_bound(_mesh->contacts->filledCells.begin(), _mesh->contacts->filledCells.end(), xsize * ysize * zmax);
		std::vector<size_t> distribution = tarray<size_t>::distribute(threads, end - begin);
		for (size_t t = 0; t <= threads; t++) {
			distribution[t] += begin - _mesh->contacts->filledCells.begin();
		}

		std::vector<std::vector<esint> > filled(threads);
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			std::vector<esint> tfilled;
			esint y, x;
			for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {
				y = _mesh->contacts->filledCells[i] % (xsize * ysize) / xsize;
				x = _mesh->contacts->filledCells[i] % xsize;
				if (xmin <= x && x < xmax && ymin <= y && y < ymax) {
					tfilled.push_back(i);
				}
			}

			filled[t].swap(tfilled);
		}

		sFilled[n].clear();
		for (size_t t = 0; t < threads; t++) {
			sFilled[n].insert(sFilled[n].end(), filled[t].begin(), filled[t].end());
		}
	}

	_mesh->contacts->gnsurface.resize(_mesh->contacts->gneighbors.size());

	for (size_t n = 0; n < _mesh->contacts->gneighbors.size(); n++) {
		std::vector<size_t> distribution = tarray<size_t>::distribute(threads, sFilled[n].size());
		sBlock[n].resize(sFilled[n].size() + 1);

		std::vector<std::vector<esint> > tIDs(threads);

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			std::vector<esint> IDs;

			if (distribution[t] < distribution[t + 1]) {
				esint prevBlock = sFilled[n][distribution[t]];
				auto cell = _mesh->contacts->grid->cbegin() + prevBlock;
				for (size_t i = distribution[t]; i < distribution[t + 1]; ++i) {
					cell += sFilled[n][i] - prevBlock;
					sFilled[n][i] = _mesh->contacts->filledCells[prevBlock = sFilled[n][i]];
					IDs.insert(IDs.end(), cell->begin(), cell->end());
					sBlock[n][i + 1] = IDs.size();
				}

				tIDs[t].swap(IDs);
			}
		}

		utils::threadDistributionToFullDistribution(sBlock[n], distribution);
		for (size_t t = 0; t < threads; t++) {
			distribution[t] = sIDs[0].size();
			sIDs[n].insert(sIDs[n].end(), tIDs[t].begin(), tIDs[t].end());
		}
		distribution[threads] = sIDs[0].size();

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			utils::sortAndRemoveDuplicity(tIDs[t]);
		}

		for (size_t t = 1; t < threads; t++) {
			tIDs[0].insert(tIDs[0].end(), tIDs[t].begin(), tIDs[t].end());
		}
		utils::sortAndRemoveDuplicity(tIDs[0]);
		_mesh->contacts->gnsurface[n].swap(tIDs[0]);

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (size_t i = distribution[t]; i < distribution[t + 1]; ++i) {
				sIDs[n][i] = std::lower_bound(_mesh->contacts->gnsurface[n].begin(), _mesh->contacts->gnsurface[n].end(), sIDs[0][i]) - _mesh->contacts->gnsurface[n].begin();
			}
		}

		distribution = tarray<size_t>::distribute(threads, _mesh->contacts->gnsurface[n].size());

		std::vector<std::vector<esint> > dist(threads);
		std::vector<std::vector<Point> > data(threads);

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			std::vector<esint> tdist;
			std::vector<Point> tdata;
			if (t == 0) {
				tdist.reserve(distribution[t + 1] - distribution[t] + 1);
				tdist.push_back(0);
			} else {
				tdist.reserve(distribution[t + 1] - distribution[t]);
			}
			if (distribution[t] < distribution[t + 1]) {
				esint prev = _mesh->contacts->gnsurface[n][distribution[t]];
				auto face = _mesh->contacts->elements->cbegin() + prev;
				for (size_t i = distribution[t]; i < distribution[t + 1]; prev = _mesh->contacts->gnsurface[n][i++]) {
					face += _mesh->contacts->gnsurface[n][i] - prev;
					tdata.insert(tdata.end(), face->begin(), face->end());
					tdist.push_back(tdata.size());
				}
			}
			dist[t].swap(tdist);
			data[t].swap(tdata);
		}
		utils::threadDistributionToFullDistribution(dist);

		for (size_t t = 0; t < threads; t++) {
			sDist[n].insert(sDist[n].end(), dist[t].begin(), dist[t].end());
			sData[n].insert(sData[n].end(), data[t].begin(), data[t].end());
		}
	}

	_mesh->contacts->gnfilled.resize(_mesh->contacts->gneighbors.size());
	_mesh->contacts->gngrid.resize(_mesh->contacts->gneighbors.size());
	_mesh->contacts->gnelements.resize(_mesh->contacts->gneighbors.size());

	if (!Communication::exchangeUnknownSize(sFilled, _mesh->contacts->gnfilled, _mesh->contacts->gneighbors)) {
		eslog::error("ESPRESO internal error: contacts - filled boxed.\n");
	}
	if (!Communication::exchangeUnknownSize(sBlock, rBlock, _mesh->contacts->gneighbors)) {
		eslog::error("ESPRESO internal error: contacts - blocks sizes.\n");
	}
	if (!Communication::exchangeUnknownSize(sIDs, rIDs, _mesh->contacts->gneighbors)) {
		eslog::error("ESPRESO internal error: contacts - IDs.\n");
	}
	if (!Communication::exchangeUnknownSize(sDist, rDist, _mesh->contacts->gneighbors)) {
		eslog::error("ESPRESO internal error: contacts - faces distribution.\n");
	}
	if (!Communication::exchangeUnknownSize(sData, rData, _mesh->contacts->gneighbors)) {
		eslog::error("ESPRESO internal error: contacts - faces data.\n");
	}

	for (size_t n = 0; n < _mesh->contacts->gneighbors.size(); ++n) {
		_mesh->contacts->gngrid[n] = new serializededata<esint, esint>(rBlock[n], rIDs[n]);
		_mesh->contacts->gnelements[n] = new serializededata<esint, Point>(rDist[n], rData[n]);
	}

	std::vector<std::vector<esint> > closeElementsDist(threads), closeElementsData(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<esint> tcloseElementsDist, tcloseElementsData;
		if (t == 0) {
			tcloseElementsDist.reserve(_mesh->contacts->elements->cbegin(t)->begin() - _mesh->contacts->elements->cbegin(t)->begin() + 1);
			tcloseElementsDist.push_back(0);
		} else {
			tcloseElementsDist.reserve(_mesh->contacts->elements->cbegin(t)->begin() - _mesh->contacts->elements->cbegin(t)->begin());
		}
		Point min, max;
		int xmin, xmax, ymin, ymax, zmin, zmax;
		size_t prevsize;
		std::vector<esint>::const_iterator zbegin, zend, ybegin, yend, xbegin, xend;
		for (auto e = _mesh->contacts->elements->cbegin(t); e != _mesh->contacts->elements->cend(t); ++e) {
			min = max = e->front();
			for (auto n = e->begin() + 1; n != e->end(); ++n) {
				min.x = std::min(min.x, n->x);
				min.y = std::min(min.y, n->y);
				min.z = std::min(min.z, n->z);
				max.x = std::max(max.x, n->x);
				max.y = std::max(max.y, n->y);
				max.z = std::max(max.z, n->z);
			}

			xmin = std::max((int)std::floor((min.x - _mesh->contacts->boundingBox[0].x - _mesh->contacts->eps - epsilon) / boxsize), 0);
			ymin = std::max((int)std::floor((min.y - _mesh->contacts->boundingBox[0].y - _mesh->contacts->eps - epsilon) / boxsize), 0);
			zmin = std::max((int)std::floor((min.z - _mesh->contacts->boundingBox[0].z - _mesh->contacts->eps - epsilon) / boxsize), 0);
			xmax = std::min((int)std::ceil ((max.x - _mesh->contacts->boundingBox[0].x + _mesh->contacts->eps + epsilon) / boxsize), xsize);
			ymax = std::min((int)std::ceil ((max.y - _mesh->contacts->boundingBox[0].y + _mesh->contacts->eps + epsilon) / boxsize), ysize);
			zmax = std::min((int)std::ceil ((max.z - _mesh->contacts->boundingBox[0].z + _mesh->contacts->eps + epsilon) / boxsize), zsize);

			auto begin = _mesh->contacts->filledCells.begin();
			auto end   = std::lower_bound(_mesh->contacts->filledCells.begin(), _mesh->contacts->filledCells.end(), xsize * ysize * (zmax - 1) + xsize * (ymax - 1) + (xmax - 1) + 1);
			auto faces = _mesh->contacts->grid->cbegin();

			auto prevcell = begin;
			prevsize = tcloseElementsData.size();
			for (int z = zmin; begin != end && z < zmax; ++z) {
				for (int y = ymin; begin != end && y < ymax; ++y) {
					begin = std::lower_bound(begin, end, xsize * ysize * z + xsize * y + xmin);
					if (begin != end) {
						faces += begin - prevcell;
						prevcell = begin;
					} else {
						break;
					}
					while (*begin < xsize * ysize * z + xsize * y + xmax) {
						tcloseElementsData.insert(tcloseElementsData.end(), faces->begin(), faces->end());
						++begin;
						if (begin != end) {
							faces += begin - prevcell;
							prevcell = begin;
						} else {
							break;
						}
					}
				}
			}
			utils::sortAndRemoveDuplicity(tcloseElementsData, prevsize);
			tcloseElementsDist.push_back(tcloseElementsData.size());
		}

		closeElementsDist[t].swap(tcloseElementsDist);
		closeElementsData[t].swap(tcloseElementsData);
	}

	utils::threadDistributionToFullDistribution(closeElementsDist);

	_mesh->contacts->closeElements = new serializededata<esint, esint>(closeElementsDist, closeElementsData);
	_mesh->contacts->gncloseElements.resize(_mesh->contacts->gneighbors.size());
	for (size_t n = 0; n < _mesh->contacts->gneighbors.size(); n++) {
		std::vector<size_t> distribution = tarray<size_t>::distribute(threads, _mesh->contacts->gnsurface[n].size());
		std::vector<std::vector<esint> > ncloseElementsDist(threads), ncloseElementsData(threads);

		double nboxsize = _mesh->contacts->eps * rGroupSize[n].front();
		int nxsize = std::ceil((boxes[2 * _mesh->contacts->gneighbors[n] + 1].x - boxes[2 * _mesh->contacts->gneighbors[n]].x) / nboxsize);
		int nysize = std::ceil((boxes[2 * _mesh->contacts->gneighbors[n] + 1].y - boxes[2 * _mesh->contacts->gneighbors[n]].y) / nboxsize);
		int nzsize = std::ceil((boxes[2 * _mesh->contacts->gneighbors[n] + 1].z - boxes[2 * _mesh->contacts->gneighbors[n]].z) / nboxsize);

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			Point min, max;
			int xmin, xmax, ymin, ymax, zmin, zmax;
			size_t prevsize;
			std::vector<esint> tcloseElementsDist, tcloseElementsData;
			if (t == 0) {
				tcloseElementsDist.reserve(distribution[t + 1] - distribution[t] + 1);
				tcloseElementsDist.push_back(0);
			} else {
				tcloseElementsDist.reserve(distribution[t + 1] - distribution[t]);
			}

			if (distribution[t] != distribution[t + 1]) {
				esint prev = _mesh->contacts->gnsurface[n][distribution[t]];
				auto e = _mesh->contacts->elements->cbegin() + prev;
				for (size_t i = distribution[t]; i < distribution[t + 1]; prev = _mesh->contacts->gnsurface[n][i++]) {
					e += _mesh->contacts->gnsurface[n][i] - prev;
					min = max = e->front();
					for (auto nn = e->begin() + 1; nn != e->end(); ++nn) {
						min.x = std::min(min.x, nn->x);
						min.y = std::min(min.y, nn->y);
						min.z = std::min(min.z, nn->z);
						max.x = std::max(max.x, nn->x);
						max.y = std::max(max.y, nn->y);
						max.z = std::max(max.z, nn->z);
					}

					xmin = std::max((int)std::floor((min.x - boxes[2 * _mesh->contacts->gneighbors[n]].x - _mesh->contacts->eps - epsilon) / nboxsize), 0);
					ymin = std::max((int)std::floor((min.y - boxes[2 * _mesh->contacts->gneighbors[n]].y - _mesh->contacts->eps - epsilon) / nboxsize), 0);
					zmin = std::max((int)std::floor((min.z - boxes[2 * _mesh->contacts->gneighbors[n]].z - _mesh->contacts->eps - epsilon) / nboxsize), 0);
					xmax = std::min((int)std::ceil ((max.x - boxes[2 * _mesh->contacts->gneighbors[n]].x + _mesh->contacts->eps + epsilon) / nboxsize), nxsize);
					ymax = std::min((int)std::ceil ((max.y - boxes[2 * _mesh->contacts->gneighbors[n]].y + _mesh->contacts->eps + epsilon) / nboxsize), nysize);
					zmax = std::min((int)std::ceil ((max.z - boxes[2 * _mesh->contacts->gneighbors[n]].z + _mesh->contacts->eps + epsilon) / nboxsize), nzsize);

					auto begin = _mesh->contacts->gnfilled[n].begin();
					auto end   = std::lower_bound(_mesh->contacts->gnfilled[n].begin(), _mesh->contacts->gnfilled[n].end(), nxsize * nysize * (zmax - 1) + nxsize * (ymax - 1) + (xmax - 1) + 1);
					auto faces = _mesh->contacts->gngrid[n]->cbegin();

					auto prevcell = begin;
					prevsize = tcloseElementsData.size();
					for (int z = zmin; begin != end && z < zmax; ++z) {
						for (int y = ymin; begin != end && y < ymax; ++y) {
							begin = std::lower_bound(begin, end, nxsize * nysize * z + nxsize * y + xmin);
							if (begin != end) {
								faces += begin - prevcell;
								prevcell = begin;
							} else {
								break;
							}
							while (*begin < nxsize * nysize * z + nxsize * y + xmax) {
								tcloseElementsData.insert(tcloseElementsData.end(), faces->begin(), faces->end());
								++begin;
								if (begin != end) {
									faces += begin - prevcell;
									prevcell = begin;
								} else {
									break;
								}
							}
						}
					}
					utils::sortAndRemoveDuplicity(tcloseElementsData, prevsize);
					tcloseElementsDist.push_back(tcloseElementsData.size());
				}
			}

			ncloseElementsDist[t].swap(tcloseElementsDist);
			ncloseElementsData[t].swap(tcloseElementsData);
		}

		utils::threadDistributionToFullDistribution(ncloseElementsDist);

		_mesh->contacts->gncloseElements[n] = new serializededata<esint, esint>(ncloseElementsDist, ncloseElementsData);
	}

	eslog::endln("MESH: CONTACT INTERFACE SEARCHED");
}



