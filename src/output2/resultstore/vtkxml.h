
#ifndef SRC_OUTPUT2_RESULTSTORE_VTKXML_H_
#define SRC_OUTPUT2_RESULTSTORE_VTKXML_H_

#include "../resultstore.h"

namespace espreso {

namespace output {

class VTKXML: public ResultStore {

public:
	virtual ~VTKXML() {};

protected:
	VTKXML(const OutputConfiguration &output, const Mesh *mesh, const std::string &path);

	virtual void store(const std::string &name, std::vector<double> &coordinates, std::vector<eslocal> &elementsTypes, std::vector<eslocal> &elementsNodes, std::vector<eslocal> &elements, DataArrays &data);
	virtual void compose(const std::string &name, const std::vector<std::string> &names);

	virtual void regionPreprocessing(const espreso::Region &region, std::vector<double> &coordinates, std::vector<eslocal> &elementsTypes, std::vector<eslocal> &elementsNodes, std::vector<eslocal> &elements);

	virtual std::string format() const =0;
	virtual void store(std::ostream &os, const std::vector<eslocal> &data) =0;
	virtual void store(std::ostream &os, const std::vector<double> &data) =0;
};

}
}


#endif /* SRC_OUTPUT2_RESULTSTORE_VTKXML_H_ */
