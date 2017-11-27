
#include "mesh.h"

#include "store/elementstore.h"
#include "store/nodestore.h"
#include "store/elementsregionstore.h"
#include "store/boundaryregionstore.h"

#include "preprocessing/meshpreprocessing.h"

#include "elements/elements.h"

#include "../basis/utilities/utils.h"
#include "../basis/utilities/communication.h"
#include "../basis/utilities/parser.h"

#include "../config/ecf/ecf.h"
#include "../config/ecf/environment.h"
#include "../config/ecf/physics/physics.h"

#include "../assembler/step.h"

#include "../old/mesh/structures/mesh.h"
#include "../old/mesh/structures/coordinates.h"
#include "../old/mesh/structures/region.h"
#include "../old/mesh/structures/elementtypes.h"
#include "../old/mesh/elements/element.h"

#include <iostream>
#include <vector>

#include "../basis/containers/serializededata.h"
#include "../basis/containers/tarray.h"
#include "../output/solution/visualization/vtklegacy.h"


using namespace espreso;


Mesh::Mesh(const ECFConfiguration &configuration)
: elements(new ElementStore(_eclasses)), nodes(new NodeStore()),
  halo(new ElementStore(_eclasses)),
  preprocessing(new MeshPreprocessing(this)),

  configuration(configuration),
  _eclasses(environment->OMP_NUM_THREADS),
  mesh(new OldMesh())
{
	size_t threads = environment->OMP_NUM_THREADS;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<Element> eclasses;
		_eclasses[t] = new Element[static_cast<int>(Element::CODE::SIZE)];

		eclasses.push_back(Point1::create());

		eclasses.push_back(Line2::create());
		eclasses.push_back(Line3::create());

		eclasses.push_back(Triangle3::create());
		eclasses.push_back(Triangle6::create());
		eclasses.push_back(Square4::create());
		eclasses.push_back(Square8::create());

		eclasses.push_back(Tetrahedron4::create(t, _eclasses[t]));
		eclasses.push_back(Tetrahedron10::create());
		eclasses.push_back(Pyramid5::create());
		eclasses.push_back(Pyramid13::create());
		eclasses.push_back(Prisma6::create());
		eclasses.push_back(Prisma15::create());
		eclasses.push_back(Hexahedron8::create(t, _eclasses[t]));
		eclasses.push_back(Hexahedron20::create());

		std::sort(eclasses.begin(), eclasses.end(), [] (const Element &e1, const Element &e2) { return static_cast<int>(e1.code) < static_cast<int>(e2.code); });

		memcpy(_eclasses[t], eclasses.data(), eclasses.size() * sizeof(Element));
	}
}

ElementsRegionStore* Mesh::eregion(const std::string &name)
{
	for (size_t r = 0; r < elementsRegions.size(); r++) {
		if (StringCompare::caseInsensitiveEq(elementsRegions[r]->name, name)) {
			return elementsRegions[r];
		}
	}
	ESINFO(ERROR) << "ESPRESO internal error: request for unknown region '" << name << "'.";
	return NULL;
}

BoundaryRegionStore* Mesh::bregion(const std::string &name)
{
	for (size_t r = 0; r < boundaryRegions.size(); r++) {
		if (StringCompare::caseInsensitiveEq(boundaryRegions[r]->name, name)) {
			return boundaryRegions[r];
		}
	}
	ESINFO(ERROR) << "ESPRESO internal error: request for unknown region '" << name << "'.";
	return NULL;
}

void Mesh::load()
{
	size_t threads = environment->OMP_NUM_THREADS;

	neighbours = mesh->neighbours();

	// LOAD NODES
	{
		std::vector<size_t> distribution = tarray<eslocal>::distribute(threads, mesh->nodes().size());
		std::vector<std::vector<Point> > coordinates(threads);
		std::vector<std::vector<eslocal> > IDs(threads);
		std::vector<std::vector<eslocal> > ranksBoundaries(threads);
		std::vector<std::vector<int> > ranksData(threads);

		ranksBoundaries.front().push_back(0);
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			size_t offset = 0;
			for (size_t n = distribution[t]; n < distribution[t + 1]; n++) {
				coordinates[t].push_back(mesh->coordinates()[n]);
				IDs[t].push_back(mesh->coordinates().globalIndex(n));
				ranksBoundaries[t].push_back(offset = offset + mesh->nodes()[n]->clusters().size());
				for (size_t c = 0; c < mesh->nodes()[n]->clusters().size(); c++) {
					ranksData[t].push_back(mesh->nodes()[n]->clusters()[c]);
				}
			}
		}

		Esutils::threadDistributionToFullDistribution(ranksBoundaries);

		nodes->size = mesh->nodes().size();
		nodes->distribution = distribution;

		nodes->IDs = new serializededata<eslocal, esglobal>(1, IDs);

		nodes->coordinates = new serializededata<eslocal, Point>(1, coordinates);
		nodes->ranks = new serializededata<eslocal, int>(ranksBoundaries, ranksData);
	}

	{
		size_t esize = mesh->elements().size();
		Communication::exscan(esize);
		std::vector<size_t> distribution = tarray<eslocal>::distribute(threads, mesh->elements().size());
		std::vector<std::vector<eslocal> > boundaries(threads), indices(threads);
		std::vector<std::vector<Element*> > epointers(threads);
		std::vector<std::vector<eslocal> > eIDs(threads);
		std::vector<std::vector<int> > body(threads), material(threads);

		boundaries.front().push_back(0);
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			size_t offset = 0;
			for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {
				switch (mesh->elements()[e]->vtkCode()) {
				case  3: epointers[t].push_back(&_eclasses[t][static_cast<int>(Element::CODE::LINE2)]); break;
				case  4: epointers[t].push_back(&_eclasses[t][static_cast<int>(Element::CODE::LINE3)]); break;

				case  5: epointers[t].push_back(&_eclasses[t][static_cast<int>(Element::CODE::TRIANGLE3)]); break;
				case  9: epointers[t].push_back(&_eclasses[t][static_cast<int>(Element::CODE::SQUARE4)]); break;
				case 22: epointers[t].push_back(&_eclasses[t][static_cast<int>(Element::CODE::TRIANGLE6)]); break;
				case 23: epointers[t].push_back(&_eclasses[t][static_cast<int>(Element::CODE::SQUARE8)]); break;

				case 10: epointers[t].push_back(&_eclasses[t][static_cast<int>(Element::CODE::TETRA4)]); break;
				case 12: epointers[t].push_back(&_eclasses[t][static_cast<int>(Element::CODE::HEXA8)]); break;
				case 13: epointers[t].push_back(&_eclasses[t][static_cast<int>(Element::CODE::PRISMA6)]); break;
				case 14: epointers[t].push_back(&_eclasses[t][static_cast<int>(Element::CODE::PYRAMID5)]); break;

				case 24: epointers[t].push_back(&_eclasses[t][static_cast<int>(Element::CODE::TETRA10)]); break;
				case 25: epointers[t].push_back(&_eclasses[t][static_cast<int>(Element::CODE::HEXA20)]); break;
				case 26: epointers[t].push_back(&_eclasses[t][static_cast<int>(Element::CODE::PRISMA15)]); break;
				case 27: epointers[t].push_back(&_eclasses[t][static_cast<int>(Element::CODE::PYRAMID13)]); break;
				}
				boundaries[t].push_back(offset = offset + mesh->elements()[e]->nodes());
				for (size_t n = 0; n < mesh->elements()[e]->nodes(); n++) {
					indices[t].push_back(mesh->elements()[e]->node(n));
				}

				eIDs[t].push_back(e + esize);
				body[t].push_back(mesh->elements()[e]->param(OldElement::Params::BODY));
				material[t].push_back(mesh->elements()[e]->param(OldElement::Params::MATERIAL));
			}
		}

		Esutils::threadDistributionToFullDistribution(boundaries);

		elements->size = mesh->elements().size();
		elements->distribution = distribution;

		elements->IDs = new serializededata<eslocal, esglobal>(1, eIDs);
		elements->nodes = new serializededata<eslocal, eslocal>(std::move(tarray<eslocal>(boundaries)), std::move(tarray<eslocal>(indices)));

		elements->body = new serializededata<eslocal, esglobal>(1, body);
		elements->material = new serializededata<eslocal, esglobal>(1, material);
		elements->epointers = new serializededata<eslocal, Element*>(1, std::move(tarray<Element*>(epointers)));
	}

	for (size_t r = 2; r < mesh->regions().size(); r++) {
		std::vector<size_t> tdistributions = tarray<size_t>::distribute(threads, mesh->regions()[r]->elements().size());
		std::vector<std::vector<eslocal> > rdistribution(threads), rdata(threads);

		if (mesh->regions()[r]->eType == ElementType::NODES) {
			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				for (size_t e = tdistributions[t]; e < tdistributions[t + 1]; e++) {
					rdata[t].push_back(mesh->regions()[r]->elements()[e]->node(0));
				}
			}
		} else {
			rdistribution[0].push_back(0);
			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				for (size_t e = tdistributions[t]; e < tdistributions[t + 1]; e++) {
					for (size_t n = 0; n < mesh->regions()[r]->elements()[e]->nodes(); n++) {
						rdata[t].push_back(mesh->regions()[r]->elements()[e]->node(n));
					}
					rdistribution[t].push_back(rdata[t].size());
				}
			}
		}

		switch (mesh->regions()[r]->eType) {
		case ElementType::ELEMENTS:
			std::cout << "region: " << mesh->regions()[r]->name << " of elements\n";
			break;
		case ElementType::FACES:
			std::cout << "region: " << mesh->regions()[r]->name << " of faces\n";
			break;
		case ElementType::EDGES:
			std::cout << "region: " << mesh->regions()[r]->name << " of edges\n";
			break;
		case ElementType::NODES:
			boundaryRegions.push_back(new BoundaryRegionStore(mesh->regions()[r]->name));
			boundaryRegions.back()->nodes = new serializededata<eslocal, eslocal>(1, rdata);
			std::sort(boundaryRegions.back()->nodes->datatarray().begin(), boundaryRegions.back()->nodes->datatarray().end());
			break;
		}
	}

	update();


//	Transformation::reclusterize(*this);
//	Transformation::partitiate(*this, 2, TFlags::SEPARATE::MATERIALS | TFlags::SEPARATE::ETYPES);
//	Transformation::computeProcessBoundaries(*this);
//	Transformation::computeDomainsBoundaries(*this); //, TFlags::ELEVEL::FACE | TFlags::ELEVEL::NODE);

//	NewOutput::VTKLegacy("nodes", _nodes, _domains);
//
//	NewOutput::VTKLegacy("processBoundaries", _processBoundaries, _nodes, false);
//	NewOutput::VTKLegacy("domainsBoundaries", _domainsBoundaries, _nodes, true);
//
//	for (size_t r = 0; r < _regions.size(); ++r) {
//		NewOutput::VTKLegacy(_regions[r]->name, _nodes, _regions[r]);
//	}
}


void Mesh::update()
{
	materials.clear();
	for (auto mat = configuration.getPhysics()->materials.begin(); mat != configuration.getPhysics()->materials.end(); ++mat) {
		materials.push_back(&mat->second);
	}

	// preprocessing->reclusterize();
	preprocessing->partitiate(3, true, true);
}

bool Mesh::prepareSolutionForOutput(const Step &step)
{
	// TODO: NUMA + load balancing

	std::vector<int> allranks(neighbours);
	allranks.push_back(environment->MPIrank);
	std::sort(allranks.begin(), allranks.end());

	auto n2i = [ & ] (int neighbour) {
		return std::lower_bound(neighbours.begin(), neighbours.end(), neighbour) - neighbours.begin();
	};

	auto r2i = [ & ] (int rank) {
		return std::lower_bound(allranks.begin(), allranks.end(), rank) - allranks.begin();
	};

	auto doffset = [&] (eslocal d, eslocal i) {
		return std::lower_bound(nodes->dintervals[d].begin(), nodes->dintervals[d].end(), i, [] (const DomainInterval &interval, eslocal i) { return interval.pindex < i; })->DOFOffset;
	};

	for (auto datait = nodes->data.begin(); datait != nodes->data.end(); ++datait) {
		NodeData* data = *datait;
		if (data->names.size()) {
			data->gathredData->clear();
			data->gathredData->resize(nodes->uniqueSize);


			std::vector<eslocal> soffsets(nodes->pintervals.size()), scouters(allranks.size());
			std::vector<std::vector<double> > sBuffer(neighbours.size()), rBuffer(neighbours.size());

			for (size_t i = 0; i < nodes->pintervals.size(); ++i) {
				soffsets[i] = scouters[r2i(nodes->pintervals[i].sourceProcess)];
				scouters[r2i(nodes->pintervals[i].sourceProcess)] += nodes->pintervals[i].end - nodes->pintervals[i].begin;
			}

			for (size_t n = 0; n < neighbours.size(); ++n) {
				if (neighbours[n] < environment->MPIrank) {
					sBuffer[n].resize(scouters[n]);
				}
			}

			#pragma omp parallel for
			for (size_t i = 0; i < nodes->pintervals.size(); ++i) {
				auto domains = nodes->idomains->cbegin() + i;
				eslocal offset, goffset;

				if (nodes->pintervals[i].sourceProcess < environment->MPIrank) {
					for (auto d = domains->begin(); d != domains->end(); ++d) {
						if (elements->firstDomain <= *d && *d < elements->firstDomain + elements->ndomains) {
							offset = doffset(*d - elements->firstDomain, i);
							for (eslocal n = nodes->pintervals[i].begin; n < nodes->pintervals[i].end; ++n, ++offset) {
								sBuffer[n2i(nodes->pintervals[i].sourceProcess)][soffsets[i]++] += (*data->decomposedData)[*d - elements->firstDomain][offset];
							}
						}
					}
				}
			}

			if (!Communication::receiveUpperUnknownSize(sBuffer, rBuffer, neighbours)) {
				ESINFO(ERROR) << "ESPRESO internal error: gather results\n";
			}

			#pragma omp parallel for
			for (size_t i = 0; i < nodes->pintervals.size(); ++i) {
				auto domains = nodes->idomains->cbegin() + i;
				eslocal offset, goffset;

				if (nodes->pintervals[i].sourceProcess == environment->MPIrank) {
					for (auto d = domains->begin(); d != domains->end(); ++d) {
						if (*d < elements->firstDomain + elements->ndomains) {
							offset = doffset(*d - elements->firstDomain, i);
							goffset = nodes->pintervals[i].globalOffset - nodes->uniqueOffset;
							for (eslocal n = nodes->pintervals[i].begin; n < nodes->pintervals[i].end; ++n, ++offset, ++goffset) {
								(*data->gathredData)[goffset] += (*data->decomposedData)[*d - elements->firstDomain][offset];
							}
						} else {
							auto neigbors = nodes->ineighborOffsets->cbegin() + i;
							for (auto neigh = neigbors->begin(); neigh != neigbors->end(); ++neigh) {
								offset = neigh->offset;
								goffset = nodes->pintervals[i].globalOffset - nodes->uniqueOffset;
								for (eslocal n = nodes->pintervals[i].begin; n < nodes->pintervals[i].end; ++n, ++offset, ++goffset) {
									(*data->gathredData)[goffset] += rBuffer[n2i(neigh->process)][offset];
								}
							}
						}
					}
					goffset = nodes->pintervals[i].globalOffset - nodes->uniqueOffset;
					for (eslocal n = nodes->pintervals[i].begin; n < nodes->pintervals[i].end; ++n, ++goffset) {
						(*data->gathredData)[goffset] /= domains->size();
					}
				}
			}
		}
	}

	return true;
}

//RegionStore* Mesh::region(const std::string &name)
//{
//	for (size_t r = 0; r < _regions.size(); r++) {
//		if (StringCompare::caseInsensitiveEq(_regions[r]->name, name)) {
//			return _regions[r];
//		}
//	}
//	ESINFO(ERROR) << "Request for unknown region '" << name << "'";
//	return NULL;
//}

