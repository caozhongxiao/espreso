
#include "newmesh.h"

#include "elements/element.h"
#include "elements/elementstore.h"

#include "../basis/utilities/utils.h"

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
: _nodes(new ElementStore()), _edges(new ElementStore()), _faces(new ElementStore()), _elems(new ElementStore())
{
	size_t threads = environment->OMP_NUM_THREADS;

	// LOAD NODES
	{
		std::vector<size_t> distribution = tarray<eslocal>::distribute(threads, mesh.nodes().size());
		std::vector<std::vector<Point> > coordinates(threads);
		std::vector<std::vector<esglobal> > IDs(threads);

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (size_t n = distribution[t]; n < distribution[t + 1]; n++) {
				coordinates[t].push_back(mesh.coordinates()[n]);
				IDs[t].push_back(mesh.coordinates().globalIndex(n));
			}
		}

		_nodes->coordinates = new serializededata<eslocal, Point>(1, coordinates);
		_nodes->IDs = new serializededata<eslocal, esglobal>(1, IDs);
	}

	auto loadElements = [&] (const std::vector<Element*> &elements) {
		std::vector<size_t> distribution = tarray<eslocal>::distribute(threads, elements.size());
		std::vector<std::vector<eslocal> > boundaries(threads), indices(threads);

		boundaries[0].push_back(0);
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			size_t offset = 0;
			for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {
				boundaries[t].push_back(offset = offset + elements[e]->nodes());
				for (size_t n = 0; n < elements[e]->nodes(); n++) {
					indices[t].push_back(elements[e]->node(n));
				}
			}
		}

		std::vector<size_t> offsets;
		for (size_t t = 0; t < threads; t++) {
			offsets.push_back(boundaries[t].size() ? boundaries[t].back() : 0);
		}
		Esutils::sizesToOffsets(offsets);

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			size_t offset = offsets[t];
			for (size_t i = 0; i < boundaries[t].size(); i++) {
				boundaries[t][i] += offset;
			}
		}

		return new serializededata<eslocal, eslocal>(tarray<eslocal>(boundaries), tarray<eslocal>(indices));
	};

	_edges->nodes = loadElements(mesh.edges());
	_faces->nodes = loadElements(mesh.faces());
	_elems->nodes = loadElements(mesh.elements());

	{
		for (size_t t = 0; t < threads; t++) {
			std::cout << "NODES [" << t << "]\n";
			auto IDs = _nodes->IDs->cbegin(t);
			auto coordinates = _nodes->coordinates->cbegin(t);
			for (; IDs != _nodes->IDs->cend(t); ++IDs, ++coordinates) {
				std::cout << IDs->front() << ": " << coordinates->front() << "\n";
			}
			std::cout << "\n";
		}
	}
	{
		for (size_t t = 0; t < threads; t++) {
			std::cout << "EDGES [" << t << "]\n";
			for (auto it = _edges->nodes->cbegin(t); it != _edges->nodes->cend(t); ++it) {
				for (size_t n = 0; n < it->size(); n++) {
					std::cout << (*it)[n] << " ";
				}
				std::cout << "\n";
			}
			std::cout << "\n";
		}
	}
	{
		for (size_t t = 0; t < threads; t++) {
			std::cout << "FACES [" << t << "]\n";
			for (auto it = _faces->nodes->cbegin(t); it != _faces->nodes->cend(t); ++it) {
				for (size_t n = 0; n < it->size(); n++) {
					std::cout << (*it)[n] << " ";
				}
				std::cout << "\n";
			}
			std::cout << "\n";
		}
	}
	{
		for (size_t t = 0; t < threads; t++) {
			std::cout << "ELEMS [" << t << "]\n";
			for (auto it = _elems->nodes->cbegin(t); it != _elems->nodes->cend(t); ++it) {
				for (size_t n = 0; n < it->size(); n++) {
					std::cout << (*it)[n] << " ";
				}
				std::cout << "\n";
			}
			std::cout << "\n";
		}
	}

	exit(0);
}

