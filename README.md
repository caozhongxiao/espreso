# ESPRESO 
#### HIGHLY PARALLEL SOLVERS FOR ENGINEERING APPLICATIONS.

ESPRESO is an ExaScale PaRallel FETI SOlver developed at Czech national supercomputing centre IT4Innovations. Main focus of the development team is to create a highly efficient parallel solver which contains several FETI based algorithms including Hybrid Total FETI method suitable for parallel machines with tens or hundreds of thousands of cores. The solver is based on highly efficient communication layer on top of pure MPI.

---
# Instalation
---

###  External Dependencies

In the current version the ESPRESO can be compiled and executed on a Linux operating system only. As of now, the library requires the [Intel MKL](https://software.intel.com/en-us/intel-mkl) and [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview>) libraries. These libraries have to be installed in the system before the ESPRESO installation process.

To compile and install the ESPRESO library a Python-based framework [Waf](https://waf.io/book/) is used. The compilation process has two phases: **configuration** and **installation**.  The former configures persistent data and checks all required headers and libraries. The latter builds the library. If something in the environment is changed, it is mandatory to re-run the configuration process.

### Building the library

The configuration process is controlled by a *build.config* script. The prefered way, how to create this script is to copy it from the install directory according to a preffered compiler. For compilation with e.g. Intel Compiler copy the appropriate script and the run configuration command:
```sh
$ cp install/build.config.icpc build.config
$ ./waf configure
```

After the configuration, the library can be installed:
```sh
$ ./waf install
```

This command builds all source files and creates the library. The type of the library can be controlled by the LIB_TYPE parameter in the configuration script.

### Set up the environment

Before running the library, environment variables need to be set. The following variables should be set according to the number of CPU cores per compute node (nCores) and number of MPI processes processed per node (PPN):

 - OMP_NUM_THREADS - should be set to nCores/PPN
 - SOLVER_NUM_THREADS - should be set to nCores/PPN
 - PAR_NUM_THREADS - should be set to nCores/PPN
 - MKL_NUM_THREADS - in the current version it should be set to 1
 - CILK_NWORKERS - in the current version it should be set to 1

It is possible to set all environment variables at once usign the following script in the ESPRESO root directory:
```sh
$ . env/threading.default ${nCores/PPN}
```


### Testing the installation
The installation of the library can be validated by the included set of benchmarks which can be executed as follows:
```sh
$ nosetest benchmarks
```
If all tests pass, the library is ready to use.

---
# ESPRESO Interface
---

The library contains a simple C API that allows to utilize FETI based solvers in third party softwares. The API has been successfully tested with an open source multiphysical simulation software [Elmer FEM](https://www.csc.fi/web/elmer)  developed by [CSC - IT Center for Science](https://www.csc.fi/).

### Call the library
Usage of the library is described by the simple example. The methods documentation can be found in the library header file.

```cpp
#include "feti4i.h"

int main() {
    // Always initialize MPI before call ESPRESO!!
	MPI_Init(&argc, &argv);

	// Calling of provider should be replaced by appropriate methods in your library !!!
	espreso::APITestESPRESODataProvider provider(&argc, &argv);

	// Solving the problem by FETI4I works in 4 steps:
	//  1. Create the stiffness matrix
	//  2. Configure the ESPRESO library
	//  3. Create an instance of a problem
	//  4. Solve the instance


	// Step 1: Create the stiffness matrix.
	//       : Fill the empty matrix by element data
	FETI4IMatrix K;
	FETI4IInt    matrixType = provider.matrixType();
	FETI4IInt    indexBase = 0;

	FETI4ICreateStiffnessMatrix(&K, matrixType, indexBase); // create matrix
	for (size_t e = 0; e < provider.elements(); e++) {
		provider.addElementMatrix(K, e); // add element data
	}


	// Step 2: Configure the ESPRESO library
	//       : Set all options to default values
	//       : Change required parameters
	FETI4IInt  iopts[FETI4I_INTEGER_OPTIONS_SIZE];
	FETI4IReal ropts[FETI4I_REAL_OPTIONS_SIZE];

	FETI4ISetDefaultIntegerOptions(iopts); // set default integer options
	FETI4ISetDefaultRealOptions(ropts); // set default real options

	/* change the default values
	iopts[FETI4IIntegerOptions::FETI4I_FETI_METHOD] = 1; // set HYBRID FETI
	iopts[FETI4IIntegerOptions::FETI4I_PRECONDITIONER] = 3; // set Dirichlet preconditioner
	*/


	// Step 3: Create an instance of a problem
	//       : Compute RHS
	//       : Fill L2G, Dirichlet, and neighbours
	FETI4IInstance            instance;
	std::vector<FETI4IReal>   rhs;
	std::vector<FETI4IInt>    dirichlet_indices;
	std::vector<FETI4IReal>   dirichlet_values;
	std::vector<FETI4IInt>    l2g;
	std::vector<FETI4IMPIInt> neighbours;

	provider.computeRHS(rhs);
	provider.fillL2G(l2g);
	provider.fillDirichlet(dirichlet_indices, dirichlet_values);
	provider.fillNeighbours(neighbours);

	FETI4ICreateInstance(
			&instance,
			K, // Stiffness matrix
			rhs.size(), rhs.data(), // RHS
			l2g.data(), // local 2 global indices mapping
			neighbours.size(), neighbours.data(), // neighbours clusters
			dirichlet_indices.size(), dirichlet_indices.data(), dirichlet_values.data(), // Dirichlet boundary condition
			iopts, ropts); // FETI4I options


	// Step 4: Solve the instance
	//       : Process the solution
	std::vector<FETI4IReal> solution(rhs.size());
	FETI4ISolve(instance, solution.size(), solution.data());


	// Finish: Destroy all data
	FETI4IDestroy(K);
	FETI4IDestroy(instance);

	MPI_Finalize();
}
```

# Licence

See the LICENSE file.