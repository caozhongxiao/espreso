
#include "mpi.h"

#include "../config/ecf/ecf.h"
#include "../mesh/mesh.h"
#include "../mesh/store/elementstore.h"
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
		std::vector<size_t> sizes(mesh.elements->ndomains);
		for (size_t d = 0; d < mesh.elements->ndomains; d++) {
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


