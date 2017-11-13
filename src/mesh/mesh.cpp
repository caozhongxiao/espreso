
#include "mesh.h"
#include "output.h"

#include "elements/elementstore.h"
#include "store/domainstore.h"
#include "store/boundarystore.h"
#include "store/regionstore.h"


#include "elements/point/point1.h"

#include "elements/line/line2.h"
#include "elements/line/line3.h"

#include "elements/plane/triangle3.h"
#include "elements/plane/triangle6.h"
#include "elements/plane/square4.h"
#include "elements/plane/square8.h"

#include "elements/volume/tetrahedron4.h"
#include "elements/volume/tetrahedron10.h"
#include "elements/volume/pyramid5.h"
#include "elements/volume/pyramid13.h"
#include "elements/volume/prisma6.h"
#include "elements/volume/prisma15.h"
#include "elements/volume/hexahedron8.h"
#include "elements/volume/hexahedron20.h"

#include "transformation/transformations.h"

#include "../basis/utilities/utils.h"
#include "../basis/utilities/communication.h"

#include "../config/ecf/environment.h"

#include "../old/mesh/structures/mesh.h"
#include "../old/mesh/structures/coordinates.h"
#include "../old/mesh/structures/region.h"
#include "../old/mesh/structures/elementtypes.h"
#include "../old/mesh/elements/element.h"

#include <iostream>
#include <vector>

#include "../basis/containers/serializededata.h"
#include "../basis/containers/tarray.h"

using namespace espreso;


Mesh::Mesh()
: _nodes(new ElementStore(_eclasses)), _edges(new ElementStore(_eclasses)), _faces(new ElementStore(_eclasses)), _elems(new ElementStore(_eclasses)), _halo(new ElementStore(_eclasses)),
  _domains(new DomainStore), _domainsBoundaries(new BoundaryStore()), _processBoundaries(new BoundaryStore()),
  _eclasses(environment->OMP_NUM_THREADS),
  mesh(NULL)
{

}

void Mesh::load()
{
	_neighbours = mesh->neighbours();
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

	// LOAD NODES
	{
		std::vector<size_t> distribution = tarray<eslocal>::distribute(threads, mesh->nodes().size());
		std::vector<std::vector<Point> > coordinates(threads);
		std::vector<std::vector<esglobal> > IDs(threads);
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

		_nodes->size = mesh->nodes().size();
		_nodes->distribution = distribution;
		_nodes->coordinates = new serializededata<eslocal, Point>(1, coordinates);
		_nodes->IDs = new serializededata<eslocal, esglobal>(1, IDs);
		_nodes->ranks = new serializededata<eslocal, int>(ranksBoundaries, ranksData);
	}

	auto loadElements = [&] (ElementStore *store, const std::vector<OldElement*> &elements) {
		std::vector<size_t> distribution = tarray<eslocal>::distribute(threads, elements.size());
		std::vector<std::vector<eslocal> > boundaries(threads), indices(threads);
		std::vector<std::vector<Element*> > epointers(threads);

		boundaries.front().push_back(0);
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			size_t offset = 0;
			for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {
				switch (elements[e]->vtkCode()) {
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
				boundaries[t].push_back(offset = offset + elements[e]->nodes());
				for (size_t n = 0; n < elements[e]->nodes(); n++) {
					indices[t].push_back(elements[e]->node(n));
				}
			}
		}

		Esutils::threadDistributionToFullDistribution(boundaries);

		store->size = elements.size();
		store->distribution = distribution;
		store->epointers = new serializededata<eslocal, Element*>(1, std::move(tarray<Element*>(epointers)));
		store->nodes = new serializededata<eslocal, eslocal>(std::move(tarray<eslocal>(boundaries)), std::move(tarray<eslocal>(indices)));
	};

	loadElements(_edges, mesh->edges());
	loadElements(_faces, mesh->faces());
	loadElements(_elems, mesh->elements());

	{
		size_t esize = mesh->elements().size();
		std::vector<size_t> distribution = tarray<eslocal>::distribute(threads, mesh->elements().size());
		Communication::exscan(esize);
		std::vector<std::vector<esglobal> > eIDs(threads);
		std::vector<std::vector<int> > body(threads), material(threads);
		for (size_t t = 0; t < threads; t++) {
			for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {
				eIDs[t].push_back(e + esize);
				body[t].push_back(mesh->elements()[e]->param(OldElement::Params::BODY));
				material[t].push_back(mesh->elements()[e]->param(OldElement::Params::MATERIAL));
			}
		}
		_elems->IDs = new serializededata<eslocal, esglobal>(1, eIDs);
		_elems->body = new serializededata<eslocal, esglobal>(1, body);
		_elems->material = new serializededata<eslocal, esglobal>(1, material);
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
			_regions.push_back(new RegionStore(mesh->regions()[r]->name, TFlags::ELEVEL::NODE));
			_regions.back()->nodes = new serializededata<eslocal, eslocal>(1, rdata);
			std::sort(_regions.back()->nodes->datatarray().begin(), _regions.back()->nodes->datatarray().end());
			break;
		}
	}

//	Transformation::reclusterize(*this);
	Transformation::partitiate(*this, 1, TFlags::SEPARATE::MATERIALS | TFlags::SEPARATE::ETYPES);
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

