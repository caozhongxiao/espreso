
#include "newmesh.h"

#include "elements/element.h"
#include "elements/elementstore.h"

#include "../basis/utilities/utils.h"

#include "../config/ecf/environment.h"

#include "../mesh/structures/mesh.h"
#include "../mesh/elements/element.h"

#include <iostream>
#include <vector>

#include "../basis/containers/serializededata.h"
#include "../basis/containers/tarray.h"

using namespace espreso;


NewMesh::NewMesh(Mesh &mesh)
: _elements(new ElementStore())
{
	std::cout << "NEW MESH init\n";

	size_t threads = environment->OMP_NUM_THREADS;
	std::vector<size_t> distribution = tarray<eslocal>::distribute(threads, mesh.elements().size());
	std::vector<std::vector<eslocal> > boundaries(threads), indices(threads);
	std::vector<size_t> offsets(threads);

	boundaries[0].push_back(0);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		size_t offset = 0;
		for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {
			boundaries[t].push_back(mesh.elements()[e]->nodes());
			offset += mesh.elements()[e]->nodes();
			for (size_t n = 0; n < mesh.elements()[e]->nodes(); n++) {
				indices[t].push_back(mesh.elements()[e]->node(n));
			}
		}
		offsets[t] = offset;
	}

	Esutils::sizesToOffsets(offsets);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		size_t offset = offsets[t];
		for (size_t i = 0; i < boundaries[t].size(); i++) {
			boundaries[t][i] = offset += boundaries[t][i];
		}
	}

	_elements->indices = new serializededata<eslocal, eslocal>(4, tarray<eslocal>(indices));
//	_elements->indices = new parrayofstructures<eslocal, eslocal>(parray<eslocal>(boundaries), parray<eslocal>(indices));

	serializededata<eslocal, eslocal> *pindices = _elements->indices;
	const serializededata<eslocal, eslocal> *constpindices = _elements->indices;

	size_t i = 0;
	for (auto it = pindices->begin(); it != pindices->end(); ++it) {
		std::cout << i++ << ": ";
		for (size_t n = 0; n != it->size(); n++) {
			std::cout << (*it)[n] << " ";
		}
		std::cout << "\n";
	}

	i = 0;
	for (auto it = constpindices->begin(); it != constpindices->end(); ++it) {
		std::cout << i++ << ": ";
		for (size_t n = 0; n != it->size(); n++) {
			std::cout << (*it)[n] << " ";
		}
		std::cout << "\n";
	}

	i = 0;
	for (size_t t = 0; t < threads; t++) {
		std::cout << t << ">>>\n";
		for (auto it = pindices->begin(t); it != pindices->end(t); ++it) {
			std::cout << i++ << ": ";
			for (size_t n = 0; n != it->size(); n++) {
				(*it)[n] = 1;
				std::cout << (*it)[n] << " ";
			}
			std::cout << "\n";
		}
	}

	i = 0;
	for (size_t t = 0; t < threads; t++) {
		std::cout << t << ">>>\n";
		for (auto it = constpindices->begin(t); it != constpindices->end(t); ++it) {
			std::cout << i++ << ": ";
			for (size_t n = 0; n != it->size(); n++) {
				std::cout << (*it)[n] << " ";
			}
			std::cout << "\n";
		}
	}

	exit(0);
}

