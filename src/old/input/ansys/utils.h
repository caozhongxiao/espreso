
#ifndef INPUT_ANSYS_UTILS_H_
#define INPUT_ANSYS_UTILS_H_

#include <cstdlib>

namespace espreso {

class OldElement;

namespace input {

class AnsysUtils {
public:
	static OldElement* createElement(eslocal *indices, eslocal n, eslocal *params);
	static OldElement* createElement(eslocal *indices, eslocal n, eslocal *params, int eType);
};


}
}




#endif /* INPUT_ANSYS_UTILS_H_ */
