
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

	TimeEval eval("SURFACE");
	eval.totalTime.startWithBarrier();
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

	size_t threads = environment->OMP_NUM_THREADS;

	_mesh->contacts->boundingBox[0] = Point(_mesh->nodes->coordinates->datatarray()[0]);
	_mesh->contacts->boundingBox[1] = Point(_mesh->nodes->coordinates->datatarray()[0]);

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

	for (size_t t = 1; t < bbox.size(); t++) {
		bbox[0].x = std::min(bbox[0].x, bbox[2 * t + 0].x);
		bbox[0].y = std::min(bbox[0].y, bbox[2 * t + 0].y);
		bbox[0].z = std::min(bbox[0].z, bbox[2 * t + 0].z);
		bbox[1].x = std::max(bbox[1].x, bbox[2 * t + 1].x);
		bbox[1].y = std::max(bbox[1].y, bbox[2 * t + 1].y);
		bbox[1].z = std::max(bbox[1].z, bbox[2 * t + 1].z);
	}

	_mesh->contacts->boundingBox[0] = bbox[0];
	_mesh->contacts->boundingBox[1] = bbox[1];

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

	eslocal boxratio = _mesh->contacts->xsize * _mesh->contacts->ysize * _mesh->contacts->zsize / _mesh->contacts->surface->elements->structures();
	_mesh->contacts->groupsize = boxratio / 20;

	event.end();
	eval.addEvent(event);

	TimeEvent event2("grid");
	event2.start();

	std::vector<std::vector<std::pair<eslocal, eslocal> > > grids(threads);

	// 2. map into grid
	///////////////////
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<std::pair<eslocal, eslocal> > grid;

		int xmin, xmax, ymin, ymax, zmin, zmax;
		size_t gsize = 0;

		auto insert = [&] (Point &min, Point &max, eslocal eindex) {
			xmin = std::floor((min.x - _mesh->contacts->boundingBox[0].x) / _mesh->contacts->eps / _mesh->contacts->groupsize);
			xmax = std::ceil((max.x - _mesh->contacts->boundingBox[0].x) / _mesh->contacts->eps / _mesh->contacts->groupsize);
			ymin = std::floor((min.y - _mesh->contacts->boundingBox[0].y) / _mesh->contacts->eps / _mesh->contacts->groupsize);
			ymax = std::ceil((max.y - _mesh->contacts->boundingBox[0].y) / _mesh->contacts->eps / _mesh->contacts->groupsize);
			zmin = std::floor((min.z - _mesh->contacts->boundingBox[0].z) / _mesh->contacts->eps / _mesh->contacts->groupsize);
			zmax = std::ceil((max.z - _mesh->contacts->boundingBox[0].z) / _mesh->contacts->eps / _mesh->contacts->groupsize);

			grid.resize(gsize + (xmax - xmin + 1) * (ymax - ymin + 1) * (zmax - zmin + 1), {0, eindex});
			for (int z = zmin; z <= zmax; ++z) {
				for (int y = ymin; y <= ymax; ++y) {
					for (int x = xmin; x <= xmax; ++x, ++gsize) {
						grid[gsize].first = _mesh->contacts->xsize * _mesh->contacts->ysize * z + _mesh->contacts->xsize * y + x;
					}
				}
			}
		};

		Point min, max;
		eslocal eindex = _mesh->contacts->surface->edistribution[t];
		for (auto e = _mesh->contacts->surface->elements->cbegin(t); e != _mesh->contacts->surface->elements->cend(t); ++e, ++eindex) {
			min = max = _mesh->nodes->coordinates->datatarray()[e->front()];
			for (auto n = e->begin() + 1; n != e->end(); ++n) {
				min.x = std::min(min.x, _mesh->nodes->coordinates->datatarray()[*n].x);
				min.y = std::min(min.y, _mesh->nodes->coordinates->datatarray()[*n].y);
				min.z = std::min(min.z, _mesh->nodes->coordinates->datatarray()[*n].z);
				max.x = std::max(max.x, _mesh->nodes->coordinates->datatarray()[*n].x);
				max.y = std::max(max.y, _mesh->nodes->coordinates->datatarray()[*n].y);
				max.z = std::max(max.z, _mesh->nodes->coordinates->datatarray()[*n].z);
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
	start("search contact interfaces");

	TimeEval eval("CONTACT");
	eval.totalTime.startWithBarrier();

	TimeEvent event3("---");
	event3.start();

	event3.end();
	eval.addEvent(event3);

	eval.totalTime.endWithBarrier();
	eval.printStatsMPI();
	MPI_Barrier(environment->MPICommunicator);
	exit(0);

	finish("search contact interfaces");
}



