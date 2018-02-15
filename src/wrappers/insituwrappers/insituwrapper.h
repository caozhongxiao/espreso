
#ifndef SRC_WRAPPERS_INSITUWRAPPERS_INSITUWRAPPER_H_
#define SRC_WRAPPERS_INSITUWRAPPERS_INSITUWRAPPER_H_

class vtkUnstructuredGrid;
class vtkCPProcessor;
class vtkCPDataDescription;

namespace espreso {

struct Step;
class Mesh;

class InSituWrapper {

public:
	InSituWrapper(const Mesh &mesh);
	~InSituWrapper();

	void update(const Step &step);

protected:
	const Mesh &_mesh;

	vtkCPProcessor *_processor;
	vtkUnstructuredGrid *_VTKGrid;
	vtkCPDataDescription *_dataDescription;
	int _timeStep;
};

}



#endif /* SRC_WRAPPERS_INSITUWRAPPERS_INSITUWRAPPER_H_ */
