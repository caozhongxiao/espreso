
#include "newmesh.h"

#include "elements/elementstore.h"

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
: _nodes(new ElementStore()), _edges(new ElementStore()), _faces(new ElementStore()), _elems(new ElementStore()), _halo(new ElementStore()),
  _eclasses(environment->OMP_NUM_THREADS)
{
	size_t threads = environment->OMP_NUM_THREADS;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		_eclasses[t].push_back(Point1::create());

		_eclasses[t].push_back(Line2::create());
		_eclasses[t].push_back(Line3::create());

		_eclasses[t].push_back(Triangle3::create());
		_eclasses[t].push_back(Triangle6::create());
		_eclasses[t].push_back(Square4::create());
		_eclasses[t].push_back(Square8::create());

		_eclasses[t].push_back(Tetrahedron4::create());
		_eclasses[t].push_back(Tetrahedron10::create());
		_eclasses[t].push_back(Pyramid5::create());
		_eclasses[t].push_back(Pyramid13::create());
		_eclasses[t].push_back(Prisma6::create());
		_eclasses[t].push_back(Prisma15::create());
		_eclasses[t].push_back(Hexahedron8::create());
		_eclasses[t].push_back(Hexahedron20::create());

		std::sort(_eclasses[t].begin(), _eclasses[t].end(), [] (const NewElement &e1, const NewElement &e2) { return static_cast<int>(e1.code) < static_cast<int>(e2.code); });
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

	_neighbours = mesh.neighbours();

	for (int rank = 0; rank < environment->MPIsize; rank++) {
		if (rank == environment->MPIrank) {
//			std::cout << "RANK: " << rank << "\n";
//			{
//				for (size_t t = 0; t < threads; t++) {
//					std::cout << "NODES [" << t << "]\n";
//					auto IDs = _nodes->IDs->cbegin(t);
//					auto coordinates = _nodes->coordinates->cbegin(t);
//					auto ranks = _nodes->ranks->cbegin(t);
//					for (; IDs != _nodes->IDs->cend(t); ++IDs, ++coordinates, ++ranks) {
//						std::cout << IDs->front() << ": " << coordinates->front() << " [";
//						for (auto rank = ranks->begin(); rank != ranks->end(); ++rank) {
//							std::cout << *rank << " ";
//						}
//						std::cout << "]\n";
//					}
//					std::cout << "\n";
//				}
//			}
//			{
//				for (size_t t = 0; t < threads; t++) {
//					std::cout << "EDGES [" << t << "]\n";
//					for (auto it = _edges->nodes->cbegin(t); it != _edges->nodes->cend(t); ++it) {
//						for (size_t n = 0; n < it->size(); n++) {
//							std::cout << (*it)[n] << " ";
//						}
//						std::cout << "\n";
//					}
//					std::cout << "\n";
//				}
//			}
//			{
//				for (size_t t = 0; t < threads; t++) {
//					std::cout << "FACES [" << t << "]\n";
//					for (auto it = _faces->nodes->cbegin(t); it != _faces->nodes->cend(t); ++it) {
//						for (size_t n = 0; n < it->size(); n++) {
//							std::cout << (*it)[n] << " ";
//						}
//						std::cout << "\n";
//					}
//					std::cout << "\n";
//				}
//			}
//			{
//				for (size_t t = 0; t < threads; t++) {
//					std::cout << "ELEMS [" << t << "]\n";
//					auto nodes = _elems->nodes->cbegin(t);
//					auto IDs = _elems->IDs->cbegin(t);
//					for (; nodes != _elems->nodes->cend(t); ++nodes, ++IDs) {
//						std::cout << IDs->front() << ": ";
//						for (size_t n = 0; n < nodes->size(); n++) {
//							std::cout << _nodes->IDs->data()[(*nodes)[n]] << " ";
//						}
//						std::cout << "\n";
//					}
//					std::cout << "\n";
//				}
//			}
		}
		MPI_Barrier(environment->MPICommunicator);
	}

	// Transformation::addLinkFromTo(*this, TFlags::ELEVEL::NODE, TFlags::ELEVEL::ELEMENT);
//	Transformation::computeDual(*this);
//	Transformation::computeDecomposedDual(*this, TFlags::SEPARATE::MATERIALS | TFlags::SEPARATE::ETYPES);
//	Transformation::computeElementCenters(*this);
	Transformation::reclusterize(*this);
	Transformation::partitiate(*this, 4, TFlags::SEPARATE::MATERIALS);

	MPI_Barrier(environment->MPICommunicator);
	exit(0);
}

