
#include "mpi.h"

#include "../config/ecf/ecf.h"
#include "../mesh/mesh.h"
#include "../mesh/store/domainstore.h"
#include "../basis/logging/logging.hpp"
#include "../input/loader.h"

// TODO: MESH
// #include "../output/datastore/espresobinaryformat.h"



using namespace espreso;

int main(int argc, char** argv)
{
	MPI_Init(&argc, &argv);
	ECFConfiguration ecf(&argc, &argv);

	size_t parts;
	std::stringstream directoryTree(ecf.decomposer.parts);
	while (directoryTree >> parts) {
		std::stringstream path;
		path << ecf.decomposer.prefix << parts * environment->MPIsize;
		// TODO: MESH
		// ESPRESOBinaryFormat::prepareDirectories(path.str(), parts);
	}

	Mesh mesh;
	Loader::load(ecf, mesh, ecf.environment.MPIrank, ecf.environment.MPIsize);
	std::stringstream decomposition(ecf.decomposer.parts);
	while (decomposition >> parts) {
		std::stringstream path;
		path << ecf.decomposer.prefix << parts * environment->MPIsize;

		// TODO: MESH
		// mesh.partitiate(parts);
		ESINFO(ALWAYS_ON_ROOT) << "Mesh partitiated to " << parts * environment->MPIsize << " parts";
		std::vector<size_t> sizes(mesh._domains->size);
		for (size_t p = 0; p < mesh._domains->size; p++) {
			// TODO: MESH
			// sizes[p] = mesh.coordinates().localSize(p);
		}
		ESINFO(ALWAYS_ON_ROOT) << "Nodes in domains: " << Info::averageValues(sizes);
		// TODO: MESH
		// ESPRESOBinaryFormat::store(mesh, path.str());
		ESINFO(ALWAYS_ON_ROOT) << "Mesh partitiated to " << parts * environment->MPIsize << " parts saved";
	}

	MPI_Finalize();
}


