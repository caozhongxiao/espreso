
#ifndef SRC_WRAPPERS_INSITUWRAPPERS_INSITUWRAPPER_H_
#define SRC_WRAPPERS_INSITUWRAPPERS_INSITUWRAPPER_H_

class vtkUnstructuredGrid;
class vtkCPProcessor;
class vtkCPDataDescription;

namespace espreso {

class Mesh;

class InSituWrapper {

public:
	InSituWrapper(const Mesh &mesh);
	~InSituWrapper();

	void update();

protected:
	const Mesh &_mesh;

	vtkCPProcessor *_processor;
	vtkUnstructuredGrid *_VTKGrid;
	vtkCPDataDescription *_dataDescription;
	int _timeStep;
};

}



#endif /* SRC_WRAPPERS_INSITUWRAPPERS_INSITUWRAPPER_H_ */
