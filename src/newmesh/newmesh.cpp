
#include "newmesh.h"
#include "output.h"

#include "elements/elementstore.h"
#include "store/domainstore.h"
#include "store/boundarystore.h"


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

#include "../mesh/structures/mesh.h"
#include "../mesh/structures/coordinates.h"
#include "../mesh/elements/element.h"

#include <iostream>
#include <vector>

#include "../basis/containers/serializededata.h"
#include "../basis/containers/tarray.h"

using namespace espreso;


NewMesh::NewMesh(Mesh &mesh)
: _nodes(new ElementStore(_eclasses)), _edges(new ElementStore(_eclasses)), _faces(new ElementStore(_eclasses)), _elems(new ElementStore(_eclasses)), _halo(new ElementStore(_eclasses)),
  _domains(new DomainStore), _domainsBoundaries(new BoundaryStore()), _processBoundaries(new BoundaryStore()),
  _neighbours(mesh.neighbours()),
  _eclasses(environment->OMP_NUM_THREADS)
{
	size_t threads = environment->OMP_NUM_THREADS;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<NewElement> eclasses;
		_eclasses[t] = new NewElement[static_cast<int>(NewElement::CODE::SIZE)];

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

		std::sort(eclasses.begin(), eclasses.end(), [] (const NewElement &e1, const NewElement &e2) { return static_cast<int>(e1.code) < static_cast<int>(e2.code); });

		memcpy(_eclasses[t], eclasses.data(), eclasses.size() * sizeof(NewElement));
	}

	// LOAD NODES
	{
		std::vector<size_t> distribution = tarray<eslocal>::distribute(threads, mesh.nodes().size());
		std::vector<std::vector<Point> > coordinates(threads);
		std::vector<std::vector<esglobal> > IDs(threads);
		std::vector<std::vector<eslocal> > ranksBoundaries(threads);
		std::vector<std::vector<int> > ranksData(threads);

		ranksBoundaries.front().push_back(0);
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			size_t offset = 0;
			for (size_t n = distribution[t]; n < distribution[t + 1]; n++) {
				coordinates[t].push_back(mesh.coordinates()[n]);
				IDs[t].push_back(mesh.coordinates().globalIndex(n));
				ranksBoundaries[t].push_back(offset = offset + mesh.nodes()[n]->clusters().size());
				for (size_t c = 0; c < mesh.nodes()[n]->clusters().size(); c++) {
					ranksData[t].push_back(mesh.nodes()[n]->clusters()[c]);
				}
			}
		}

		Esutils::threadDistributionToFullDistribution(ranksBoundaries);

		_nodes->size = mesh.nodes().size();
		_nodes->distribution = distribution;
		_nodes->coordinates = new serializededata<eslocal, Point>(1, coordinates);
		_nodes->IDs = new serializededata<eslocal, esglobal>(1, IDs);
		_nodes->ranks = new serializededata<eslocal, int>(ranksBoundaries, ranksData);
	}

	auto loadElements = [&] (ElementStore *store, const std::vector<Element*> &elements) {
		std::vector<size_t> distribution = tarray<eslocal>::distribute(threads, elements.size());
		std::vector<std::vector<eslocal> > boundaries(threads), indices(threads);
		std::vector<std::vector<NewElement*> > epointers(threads);

		boundaries.front().push_back(0);
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			size_t offset = 0;
			for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {
				switch (elements[e]->vtkCode()) {
				case  3: epointers[t].push_back(&_eclasses[t][static_cast<int>(NewElement::CODE::LINE2)]); break;
				case  4: epointers[t].push_back(&_eclasses[t][static_cast<int>(NewElement::CODE::LINE3)]); break;

				case  5: epointers[t].push_back(&_eclasses[t][static_cast<int>(NewElement::CODE::TRIANGLE3)]); break;
				case  9: epointers[t].push_back(&_eclasses[t][static_cast<int>(NewElement::CODE::SQUARE4)]); break;
				case 22: epointers[t].push_back(&_eclasses[t][static_cast<int>(NewElement::CODE::TRIANGLE6)]); break;
				case 23: epointers[t].push_back(&_eclasses[t][static_cast<int>(NewElement::CODE::SQUARE8)]); break;

				case 10: epointers[t].push_back(&_eclasses[t][static_cast<int>(NewElement::CODE::TETRA4)]); break;
				case 12: epointers[t].push_back(&_eclasses[t][static_cast<int>(NewElement::CODE::HEXA8)]); break;
				case 13: epointers[t].push_back(&_eclasses[t][static_cast<int>(NewElement::CODE::PRISMA6)]); break;
				case 14: epointers[t].push_back(&_eclasses[t][static_cast<int>(NewElement::CODE::PYRAMID5)]); break;

				case 24: epointers[t].push_back(&_eclasses[t][static_cast<int>(NewElement::CODE::TETRA10)]); break;
				case 25: epointers[t].push_back(&_eclasses[t][static_cast<int>(NewElement::CODE::HEXA20)]); break;
				case 26: epointers[t].push_back(&_eclasses[t][static_cast<int>(NewElement::CODE::PRISMA15)]); break;
				case 27: epointers[t].push_back(&_eclasses[t][static_cast<int>(NewElement::CODE::PYRAMID13)]); break;
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
		store->epointers = new serializededata<eslocal, NewElement*>(1, std::move(tarray<NewElement*>(epointers)));
		store->nodes = new serializededata<eslocal, eslocal>(std::move(tarray<eslocal>(boundaries)), std::move(tarray<eslocal>(indices)));
	};

	loadElements(_edges, mesh.edges());
	loadElements(_faces, mesh.faces());
	loadElements(_elems, mesh.elements());

	{
		size_t esize = mesh.elements().size();
		std::vector<size_t> distribution = tarray<eslocal>::distribute(threads, mesh.elements().size());
		Communication::exscan(esize);
		std::vector<std::vector<esglobal> > eIDs(threads);
		std::vector<std::vector<int> > body(threads), material(threads);
		for (size_t t = 0; t < threads; t++) {
			for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {
				eIDs[t].push_back(e + esize);
				body[t].push_back(mesh.elements()[e]->param(Element::Params::BODY));
				material[t].push_back(mesh.elements()[e]->param(Element::Params::MATERIAL));
			}
		}
		_elems->IDs = new serializededata<eslocal, esglobal>(1, eIDs);
		_elems->body = new serializededata<eslocal, esglobal>(1, body);
		_elems->material = new serializededata<eslocal, esglobal>(1, material);
	}

	// Transformation::addLinkFromTo(*this, TFlags::ELEVEL::NODE, TFlags::ELEVEL::ELEMENT);
//	Transformation::computeDual(*this);
//	Transformation::computeDecomposedDual(*this, TFlags::SEPARATE::MATERIALS | TFlags::SEPARATE::ETYPES);
//	Transformation::computeElementCenters(*this);
	Transformation::reclusterize(*this);
	Transformation::partitiate(*this, 4, TFlags::SEPARATE::MATERIALS | TFlags::SEPARATE::ETYPES);
	Transformation::computeDomainsBoundaries(*this); //, TFlags::ELEVEL::FACE | TFlags::ELEVEL::NODE);

	NewOutput::VTKLegacy("test", _domainsBoundaries, _nodes);

	MPI_Barrier(environment->MPICommunicator);
	MPI_Finalize();
	exit(0);
}

