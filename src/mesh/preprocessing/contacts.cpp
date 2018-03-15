
#include "meshpreprocessing.h"

#include "../mesh.h"
#include "../elements/element.h"

#include "../store/elementstore.h"
#include "../store/nodestore.h"
#include "../store/surfacestore.h"
#include "../store/contactstore.h"

#include "../../basis/containers/serializededata.h"
#include "../../basis/utilities/communication.h"
#include "../../basis/utilities/utils.h"
#include "../../basis/logging/timeeval.h"

#include "../../config/ecf/environment.h"

#include <algorithm>
#include <numeric>

using namespace espreso;

void MeshPreprocessing::computeBodiesSurface()
{
	start("computation surface");

	if (_mesh->elements->neighbors == NULL) {
		this->computeElementsNeighbors();
	}

	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<std::vector<eslocal> > faces(threads), facesDistribution(threads), ecounters(threads, std::vector<eslocal>((int)Element::CODE::SIZE));
	std::vector<std::vector<Element*> > fpointers(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto nodes = _mesh->elements->nodes->cbegin(t);
		auto neighs = _mesh->elements->neighbors->cbegin(t);
		const auto &epointers = _mesh->elements->epointers->datatarray();

		std::vector<eslocal> fdist, fdata, ecounter((int)Element::CODE::SIZE);
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

	serializededata<eslocal, Element*>::balance(1, fpointers);
	_mesh->surface->epointers = new serializededata<eslocal, Element*>(1, fpointers);
	_mesh->surface->ecounters = ecounters[0];

	_mesh->surface->edistribution = _mesh->surface->epointers->datatarray().distribution();

	if (_mesh->surface->ecounters[(int)Element::CODE::TRIANGLE3] == _mesh->surface->edistribution.back()) {
		serializededata<eslocal, eslocal>::balance(3, faces, &_mesh->surface->edistribution);
		_mesh->surface->elements = new serializededata<eslocal, eslocal>(3, faces);
		_mesh->surface->triangles = _mesh->surface->elements;
		_mesh->surface->tdistribution = _mesh->surface->edistribution;
	} else {
		Esutils::threadDistributionToFullDistribution(facesDistribution);
		serializededata<eslocal, eslocal>::balance(facesDistribution, faces, &_mesh->surface->edistribution);
		_mesh->surface->elements = new serializededata<eslocal, eslocal>(facesDistribution, faces);
	}

	finish("computation surface");
}

void MeshPreprocessing::computeSurfaceLocations()
{
	start("computation surface location");

	size_t threads = environment->OMP_NUM_THREADS;

	TimeEval eval("SURFACE");
	eval.totalTime.startWithBarrier();

	TimeEvent event0("elements with coordinates");
	event0.start();

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
		_mesh->contacts->elements = new serializededata<eslocal, Point>(tarray<eslocal>(_mesh->surface->elements->boundarytarray()), points);
	} else {
		_mesh->contacts->elements = new serializededata<eslocal, Point>(3, points);
	}

	event0.end();
	eval.addEvent(event0);

	TimeEvent event("bounding box");
	event.start();

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
		eslocal eint = _mesh->nodes->externalIntervals[t];
		Point min = _mesh->nodes->coordinates->datatarray().front();
		Point max = _mesh->nodes->coordinates->datatarray().front();
		for (eslocal n = _mesh->nodes->pintervals[eint].begin; n < _mesh->nodes->pintervals[eint].end; ++n) {
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

	_mesh->contacts->xsize = std::round((_mesh->contacts->boundingBox[1].x - _mesh->contacts->boundingBox[0].x) / _mesh->contacts->eps);
	_mesh->contacts->ysize = std::round((_mesh->contacts->boundingBox[1].y - _mesh->contacts->boundingBox[0].y) / _mesh->contacts->eps);
	_mesh->contacts->zsize = std::round((_mesh->contacts->boundingBox[1].z - _mesh->contacts->boundingBox[0].z) / _mesh->contacts->eps);

	double avgsize = (_mesh->contacts->xsize + _mesh->contacts->ysize + _mesh->contacts->zsize) / 3.;
	_mesh->contacts->groupsize = std::max(avgsize / std::pow(_mesh->contacts->surface->elements->structures(), 1. / 3) / 1.5, 1.);

	MPI_Allreduce(&_mesh->contacts->boundingBox[0], &_mesh->contacts->globalBox[0], 3, MPI_DOUBLE, MPI_MIN, environment->MPICommunicator);
	MPI_Allreduce(&_mesh->contacts->boundingBox[1], &_mesh->contacts->globalBox[1], 3, MPI_DOUBLE, MPI_MAX, environment->MPICommunicator);

	_mesh->contacts->xbegin = std::floor((_mesh->contacts->boundingBox[0].x - _mesh->contacts->globalBox[0].x) / _mesh->contacts->eps / _mesh->contacts->groupsize);
	_mesh->contacts->xend = std::ceil((_mesh->contacts->boundingBox[1].x - _mesh->contacts->globalBox[0].x) / _mesh->contacts->eps / _mesh->contacts->groupsize);
	_mesh->contacts->ybegin = std::floor((_mesh->contacts->boundingBox[0].y - _mesh->contacts->globalBox[0].y) / _mesh->contacts->eps / _mesh->contacts->groupsize);
	_mesh->contacts->yend = std::ceil((_mesh->contacts->boundingBox[1].y - _mesh->contacts->globalBox[0].y) / _mesh->contacts->eps / _mesh->contacts->groupsize);
	_mesh->contacts->zbegin = std::floor((_mesh->contacts->boundingBox[0].z - _mesh->contacts->globalBox[0].z) / _mesh->contacts->eps / _mesh->contacts->groupsize);
	_mesh->contacts->zend = std::ceil((_mesh->contacts->boundingBox[1].z - _mesh->contacts->globalBox[0].z) / _mesh->contacts->eps / _mesh->contacts->groupsize);

	event.end();
	eval.addEvent(event);

	TimeEvent event2("grid");
	event2.start();

	std::vector<std::vector<std::pair<eslocal, eslocal> > > grids(threads);

	double boxsize = _mesh->contacts->eps * _mesh->contacts->groupsize;
	size_t xsize = _mesh->contacts->xend - _mesh->contacts->xbegin;
	size_t ysize = _mesh->contacts->yend - _mesh->contacts->ybegin;
	size_t zsize = _mesh->contacts->zend - _mesh->contacts->zbegin;

	// 2. map into grid
	///////////////////
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<std::pair<eslocal, eslocal> > grid;

		int xmin, xmax, ymin, ymax, zmin, zmax;
		size_t gsize = 0;

		auto insert = [&] (Point &min, Point &max, eslocal eindex) {
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
		eslocal eindex = _mesh->contacts->surface->edistribution[t];
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
		Esutils::sortWithInplaceMerge(grids[0], gdistribution);
	}

	event2.end();
	eval.addEvent(event2);

	// 3. create sparse structure
	///////////////////

	TimeEvent event3("sparse structure");
	event3.start();

	std::vector<size_t> distribution = tarray<std::pair<eslocal, eslocal> >::distribute(threads, grids[0].size());

	for (size_t t = 1; t < threads; t++) {
		while (distribution[t] != distribution[t + 1] && grids[0][distribution[t]].first == grids[0][distribution[t] - 1].first) {
			++distribution[t];
		}
	}

	std::vector<std::vector<eslocal> > gdist(threads), gdata(threads), gfilled(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<eslocal> dist, data, filled;
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

	Esutils::threadDistributionToFullDistribution(gdist);
	for (size_t t = 1; t < threads; t++) {
		gdist[0].insert(gdist[0].end(), gdist[t].begin(), gdist[t].end());
		gdata[0].insert(gdata[0].end(), gdata[t].begin(), gdata[t].end());
		gfilled[0].insert(gfilled[0].end(), gfilled[t].begin(), gfilled[t].end());
	}

	// never threaded
	gdist.resize(1);
	gdata.resize(1);

	_mesh->contacts->filledCells.swap(gfilled[0]);
	_mesh->contacts->grid = new serializededata<eslocal, eslocal>(gdist, gdata);

	event3.end();
	eval.addEvent(event3);

	eval.totalTime.endWithBarrier();
	eval.printStatsMPI();

	finish("computation surface location");
}

void MeshPreprocessing::searchContactInterfaces()
{
	return;
	start("search contact interfaces");

	double epsilon = 1e-6;

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

	TimeEval eval("CONTACT");
	eval.totalTime.startWithBarrier();

	TimeEvent event3("---");
	event3.start();

	std::vector<Point> boxes(2 * environment->MPIsize);
	std::vector<Point> intersections;
	MPI_Allgather(_mesh->contacts->boundingBox, 6, MPI_DOUBLE, boxes.data(), 6, MPI_DOUBLE, environment->MPICommunicator);

	for (int r = 0; r < environment->MPIsize; r++) {
		if (r != environment->MPIrank) {
			if (areIntersected(_mesh->contacts->boundingBox, boxes.data() + 2 * r)) {
				_mesh->contacts->neighbors.push_back(r);
			}
		}
	}

	std::vector<std::vector<size_t> > sGroupSize(_mesh->contacts->neighbors.size(), { _mesh->contacts->groupsize }), rGroupSize(_mesh->contacts->neighbors.size());
	std::vector<std::vector<eslocal> > sFilled(_mesh->contacts->neighbors.size()), rFilled(_mesh->contacts->neighbors.size());
	std::vector<std::vector<eslocal> > sDist(_mesh->contacts->neighbors.size()), rDist(_mesh->contacts->neighbors.size());
	std::vector<std::vector<eslocal> > sData(_mesh->contacts->neighbors.size()), rData(_mesh->contacts->neighbors.size());

	intersections.resize(2 * _mesh->contacts->neighbors.size());
	double boxsize = _mesh->contacts->eps * _mesh->contacts->groupsize;
	size_t xsize = _mesh->contacts->xend - _mesh->contacts->xbegin;
	size_t ysize = _mesh->contacts->yend - _mesh->contacts->ybegin;
	size_t zsize = _mesh->contacts->zend - _mesh->contacts->zbegin;

	#pragma omp parallel for
	for (size_t n = 0; n < _mesh->contacts->neighbors.size(); n++) {
		intersect(intersections.data() + 2 * n, _mesh->contacts->boundingBox, boxes.data() + 2 * _mesh->contacts->neighbors[n]);

		int xmin = std::floor((intersections[2 * n].x - _mesh->contacts->boundingBox[0].x - epsilon) / boxsize);
		int xmax = std::ceil((intersections[2 * n + 1].x - _mesh->contacts->boundingBox[0].x + epsilon) / boxsize);
		int ymin = std::floor((intersections[2 * n].y - _mesh->contacts->boundingBox[0].y - epsilon) / boxsize);
		int ymax = std::ceil((intersections[2 * n + 1].y - _mesh->contacts->boundingBox[0].y + epsilon) / boxsize);
		int zmin = std::floor((intersections[2 * n].z - _mesh->contacts->boundingBox[0].z - epsilon) / boxsize);
		int zmax = std::ceil((intersections[2 * n + 1].z - _mesh->contacts->boundingBox[0].z + epsilon) / boxsize);

		sDist[n].push_back(0);
		auto filled = _mesh->contacts->filledCells.begin();
		eslocal position;
		auto cell = _mesh->contacts->grid->cbegin();
		for (int z = zmin; z < zmax; ++z) {
			for (int y = ymin; y < ymax; ++y) {
				for (int x = xmin; x < xmax; ++x) {
					position = xsize * ysize * z + xsize * y + x;
					if (position < _mesh->contacts->filledCells.back()) {
						while (*filled < position) {
							++filled;
							++cell;
						}
						if (*filled == position) {
							sFilled[n].push_back(position);
							sData[n].insert(sData[n].end(), cell->begin(), cell->end());
							sDist[n].push_back(sDist[n].back() + cell->size());
						}
					}
				}
			}
		}

	}




	if (!Communication::exchangeKnownSize(sGroupSize, rGroupSize, _mesh->contacts->neighbors)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange group size.";
	}

	if (!Communication::exchangeKnownSize(sFilled, rFilled, _mesh->contacts->neighbors)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange group size.";
	}

	if (!Communication::exchangeKnownSize(sDist, rDist, _mesh->contacts->neighbors)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange group size.";
	}

	if (!Communication::exchangeKnownSize(sData, rData, _mesh->contacts->neighbors)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange group size.";
	}

	event3.end();
	eval.addEvent(event3);
	eval.totalTime.endWithBarrier();
	eval.printStatsMPI();
//	MPI_Barrier(environment->MPICommunicator);
//	exit(0);

	finish("search contact interfaces");
}



