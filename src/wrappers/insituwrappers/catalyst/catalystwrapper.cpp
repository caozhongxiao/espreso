
#include "../insituwrapper.h"

#include "../../../basis/containers/point.h"
#include "../../../basis/containers/serializededata.h"
#include "../../../basis/logging/logging.h"

#include "../../../assembler/step.h"
#include "../../../mesh/mesh.h"
#include "../../../mesh/store/elementstore.h"
#include "../../../mesh/store/nodestore.h"

#include "../../../output/result/visualization/vtkwritter.h"

#include "vtkNew.h"

#include "vtkDoubleArray.h"

#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkSmartPointer.h"
#include "vtkUnstructuredGrid.h"

#include "vtkCPProcessor.h"
#include "vtkCPDataDescription.h"
#include "vtkCPPythonScriptPipeline.h"
#include "vtkCPInputDataDescription.h"
#include "vtkFieldData.h"

using namespace espreso;

InSituWrapper::InSituWrapper(const Mesh &mesh)
: _mesh(mesh), _timeStep(0)
{
	_processor = vtkCPProcessor::New();
	_processor->Initialize();

	vtkSmartPointer<vtkCPPythonScriptPipeline> pipeline = vtkSmartPointer<vtkCPPythonScriptPipeline>::New();
	pipeline->Initialize("catalyst.py");
	_processor->AddPipeline(pipeline.GetPointer());

	_dataDescription = vtkCPDataDescription::New();
	_dataDescription->AddInput("input");
	_dataDescription->SetTimeData(0.0, 0);

	_VTKGrid = vtkUnstructuredGrid::New();

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	points->SetDataTypeToDouble();
	points->SetNumberOfPoints(_mesh.nodes->coordinates->datatarray().size());
	for (size_t i = 0; i < _mesh.nodes->coordinates->datatarray().size(); ++i) {
		points->InsertNextPoint(_mesh.nodes->coordinates->datatarray()[i].x, _mesh.nodes->coordinates->datatarray()[i].y, _mesh.nodes->coordinates->datatarray()[i].z);
	}
	_VTKGrid->SetPoints(points);

	_VTKGrid->Allocate(_mesh.elements->size);
	std::vector<vtkIdType> nodes(20);
	auto epointers = _mesh.elements->epointers->datatarray().cbegin();
	for (auto e = _mesh.elements->nodes->cbegin(); e != _mesh.elements->nodes->cend(); ++e, ++epointers) {
		nodes.clear();
		nodes.insert(nodes.end(), e->begin(), e->end());
		_VTKGrid->InsertNextCell(VTKWritter::ecode((*epointers)->code), nodes.size(), nodes.data());
	}

	_dataDescription->GetInputDescriptionByName("input")->SetGrid(_VTKGrid);
	_processor->CoProcess(_dataDescription);
}

void InSituWrapper::update(const Step &step)
{
	_dataDescription->SetTimeData(step.currentTime, ++_timeStep);

	std::vector<double> data;
	for (size_t di = 0; di < _mesh.nodes->data.size(); di++) {
		data.resize(_mesh.nodes->data[di]->dimension * _mesh.nodes->size, 1);
		if (_mesh.nodes->data[di]->names.size()) {
			vtkNew<vtkDoubleArray> vtkArray;
			vtkArray->SetName(_mesh.nodes->data[di]->names.front().c_str());
			vtkArray->SetNumberOfComponents(_mesh.nodes->data[di]->dimension);
			vtkArray->SetArray(data.data(), _mesh.nodes->size, 1);

			_VTKGrid->GetPointData()->AddArray(vtkArray.GetPointer());
			if (_VTKGrid->GetPointData()->GetNumberOfArrays() == 1) {
				_VTKGrid->GetPointData()->SetActiveScalars(_mesh.nodes->data[di]->names.front().c_str());
			}
		}
	}
	_processor->CoProcess(_dataDescription);
}

InSituWrapper::~InSituWrapper()
{
	_dataDescription->Delete();
	_processor->Finalize();
	_processor->Delete();
	_VTKGrid->Delete();
}



