#include <iostream>
#include <vector>
#include "Adaptor.h"

#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkCPDataDescription.h>
#include <vtkCPInputDataDescription.h>
#include <vtkCPProcessor.h>
#include <vtkCPPythonScriptPipeline.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>

#include<mpi.h>

#define MAX_NUMBER_OF_NODES 20

namespace
{
  vtkCPProcessor* Processor = NULL;
  vtkUnstructuredGrid* VTKGrid;

  void BuildVTKGrid(mesh::Mesh *mesh, 
										std::vector<std::vector<eslocal> > &l2g_vec,
										std::vector<double>& grid_points)
  {

	//UNDER CONSTRUCTION
	int MPIsize = 1;
	int MPIrank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);
	MPI_Comm_size(MPI_COMM_WORLD, &MPIsize);

	MPI_Barrier(MPI_COMM_WORLD);
	std::cout << "##############################################here\n";
	std::cout << "MPIsize = " << MPIsize << "\n";
	std::cout << "MPIrank = " << MPIrank << "\n";
	std::cout << "##############################################here\n";
	MPI_Barrier(MPI_COMM_WORLD);

	const std::vector<mesh::Element*> &elements = mesh->getElements();
  size_t numCells = 0;
  for(size_t iEl=0; iEl<elements.size(); iEl++) { 
		numCells+=elements[iEl]->size(); 
	}


// OLD CODE NODES    -- start ----------------------------------------------------------|
    vtkNew<vtkDoubleArray> pointArray;
    pointArray->SetNumberOfComponents(3);

    size_t numpoints = grid_points.size();
    double* array_grid_points = &grid_points[0];  //Change to array instead of vector

    pointArray->SetArray(array_grid_points, numpoints, 1);

    vtkNew<vtkPoints> points;
    points->SetData(pointArray.GetPointer());  //Need to convert vector to vtk datatype
    VTKGrid->SetPoints(points.GetPointer());
// OLD CODE NODES    -- end ------------------------------------------------------------|


// NEW CODE NODES    -- start ----------------------------------------------------------|

		size_t nSubClst = l2g_vec.size();
		size_t cnt = 0;

		size_t n_points = 0;
		for (size_t d = 0; d < l2g_vec.size(); d++) {
			n_points += l2g_vec[d].size();
		}

printf("n_points = %d \n", n_points);     // n_points is number of nodes 
printf("numpoints = %d \n", numpoints);   // numpoints is number of coordinates ( nodes * 3 )

// TODO Why it does not work ??? // 
mesh::Coordinates &coords  = mesh->coordinates();
printf(" one element of coord: %3.3e \n",coords[0].x);
/*
		size_t *coord = new size_t [n_points*3];

		for (size_t d = 0; d < nSubClst; d++) {
			mesh::Point center;
			for (size_t c = 0; c < l2g_vec[d].size(); c++) {
				center += _coordinates[l2g_vec[d][c]];
			}
			center /= l2g_vec[d].size();

			for (size_t i = 0; i < l2g_vec[d].size(); i++) {
				mesh::Point x = _coordinates[l2g_vec[d][i]];
				x = center + (x - center) * shrinking;
				coord[3*i+0] = x.x;
				coord[3*i+1] = x.y;
				coord[3*i+2] = x.z;
				}
			}
    vtkNew<vtkDoubleArray> pointArray;
    pointArray->SetNumberOfComponents(3);

    size_t numpoints = n_points*3;

    pointArray->SetArray(grid_points, numpoints, 1);

    vtkNew<vtkPoints> points;
    points->SetData(pointArray.GetPointer());  //Need to convert vector to vtk datatype
    VTKGrid->SetPoints(points.GetPointer());

// NEW CODE NODES    -- end ------------------------------------------------------------|
*/


// NEW CODE ELEMENTS -- start ----------------------------------------------------------|
    VTKGrid->Allocate(static_cast<vtkIdType>(numCells));
	  vtkIdType tmp[MAX_NUMBER_OF_NODES]; /* max number of nodes predicted to be '20' */
    for(size_t iEl=0; iEl<elements.size(); iEl++) { 
			for(size_t jNod=0; jNod<elements[iEl]->size(); jNod++){
				tmp[jNod] = elements[iEl]->node(jNod);
			}
			VTKGrid->InsertNextCell(elements[iEl]->vtkCode(), elements[iEl]->size(), &tmp[0]);
		}
// NEW CODE ELEMENTS -- end ------------------------------------------------------------|

  }


  void UpdateVTKAttributes(std::vector<std::vector<double> >& prim_solution, std::vector<float>& decomposition_values)
  {
	 int MPIsize;
	 int MPIrank;
	 MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);
	 MPI_Comm_size(MPI_COMM_WORLD, &MPIsize);

	unsigned int mycounter = 0;
//	if (MPIrank == 0 )
//	{
		for (size_t i = 0; i < prim_solution.size(); i++)
		{
			mycounter += prim_solution[i].size();
		}
//		MPI_Bcast(&mycounter, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
//	}
//	if (MPIrank != 0)
//	{
//		MPI_Bcast(&mycounter, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
//	}
//	MPI_Barrier(MPI_COMM_WORLD);

	double displacement_array[mycounter];
	unsigned int size = mycounter;

    mycounter = 0;
	for (size_t i = 0; i < prim_solution.size(); i++) {	//note prim_solution.size = number of subdomains. all elements in prim_solution[0]
		for (size_t j = 0; j < (prim_solution[i].size()/3); j++) {
			displacement_array[3*mycounter + 0] = prim_solution[i][3 * j + 0];
			displacement_array[3*mycounter + 1] = prim_solution[i][3 * j + 1];
			displacement_array[3*mycounter + 2] = prim_solution[i][3 * j + 2];

			mycounter ++;
			//std::cout << "counter: " << mycounter << " From Process " << MPIrank << "\n";
		}
	}

	 if(VTKGrid->GetPointData()->GetNumberOfArrays() == 0)
      {
		  // displacement array
		  vtkNew<vtkDoubleArray> displacement;
		  displacement->SetName("displacement");
		  displacement->SetNumberOfComponents(3);
		  //displacement->SetNumberOfTuples(static_cast<vtkIdType>(prim_solution[0].size()/3));
		  displacement->SetNumberOfTuples(static_cast<vtkIdType>(mycounter));
		  VTKGrid->GetPointData()->AddArray(displacement.GetPointer());
      }
	  vtkDoubleArray* displacement = vtkDoubleArray::SafeDownCast(
	  VTKGrid->GetPointData()->GetArray("displacement"));
	  double* displacementData = displacement_array;
	  vtkIdType numTuples = displacement->GetNumberOfTuples();
	  int counter = 0;
	  for(vtkIdType i=0;i<numTuples;i++)
	  {
	      double values[3] = {displacementData[i*3], displacementData[i*3 + 1], displacementData[i*3 + 2]};
	      displacement->SetTupleValue(counter, values);

	      counter++;
	   }

	//TODO Why decomposition_array is rebuild?
	 //DECOMPOSITION
	float* decomposition_array = &decomposition_values[0];
	float numCells = decomposition_values.size();
    if(VTKGrid->GetCellData()->GetNumberOfArrays() == 0)
      {
       //decomposition array
      vtkNew<vtkFloatArray> decomposition;
      decomposition->SetName("decomposition");
      decomposition->SetNumberOfComponents(1);
      VTKGrid->GetCellData()->AddArray(decomposition.GetPointer());
      }

    vtkFloatArray* decomposition = vtkFloatArray::SafeDownCast(
      VTKGrid->GetCellData()->GetArray("decomposition"));
    // The decomposition array is a scalar array so we can reuse
    // memory as long as we ordered the points properly.
    float* decompositionData = decomposition_array;
    decomposition->SetArray(decompositionData, static_cast<vtkIdType>(numCells), 1);
  }



  void BuildVTKDataStructures(mesh::Mesh *mesh, 
															std::vector<std::vector<eslocal> > &l2g_vec,
															std::vector<double>& grid_points, 
															 std::vector<std::vector<double> >& prim_solution, std::vector<float>& decomposition_values)
  {
    if(VTKGrid == NULL)
      {
      // The grid structure isn't changing so we only build it
      // the first time it's needed. If we needed the memory
      // we could delete it and rebuild as necessary.
      VTKGrid = vtkUnstructuredGrid::New();
      BuildVTKGrid(mesh, l2g_vec,grid_points);
      }
    UpdateVTKAttributes(prim_solution, decomposition_values);
  }
}

namespace Adaptor
{

  void Initialize(int numScripts, char* scripts[])
  {
    if(Processor == NULL)
      {
      Processor = vtkCPProcessor::New();
      Processor->Initialize();
      }
    else
      {
      Processor->RemoveAllPipelines();
      }

    //for(int i=1;i<numScripts;i++)
    //  {
      vtkNew<vtkCPPythonScriptPipeline> pipeline;
      pipeline->Initialize(scripts[numScripts-1]);
      Processor->AddPipeline(pipeline.GetPointer());
    //  }
  }

  void Finalize()
  {
    if(Processor)
      {
      Processor->Delete();
      Processor = NULL;
      }
    if(VTKGrid)
      {
      VTKGrid->Delete();
      VTKGrid = NULL;
      }
  }

  void CoProcess(	mesh::Mesh *mesh, 
									std::vector<std::vector<eslocal> > &l2g_vec,
									std::vector<double>& grid_points, 
									std::vector<std::vector<double> >& prim_solution, 
									std::vector<float>& decomposition_values, double time,
                 	unsigned int timeStep, bool lastTimeStep)
  {
    vtkNew<vtkCPDataDescription> dataDescription;
    dataDescription->AddInput("input");
    dataDescription->SetTimeData(time, timeStep);
    if(lastTimeStep == true)
      {
      // assume that we want to all the pipelines to execute if it
      // is the last time step.
      dataDescription->ForceOutputOn();
      }
    if(Processor->RequestDataDescription(dataDescription.GetPointer()) != 0)
      {
      BuildVTKDataStructures(mesh,l2g_vec, grid_points, 
														prim_solution, decomposition_values);
      dataDescription->GetInputDescriptionByName("input")->SetGrid(VTKGrid);
      Processor->CoProcess(dataDescription.GetPointer());
      }
  }
} // end of Catalyst namespace
