
#ifndef SRC_OUTPUT_RESULT_VISUALIZATION_VTKWRITTER_H_
#define SRC_OUTPUT_RESULT_VISUALIZATION_VTKWRITTER_H_

#include "../../../mesh/elements/element.h"

namespace espreso {

struct VTKWritter {

	static int ecode(Element::CODE &code)
	{
		switch (code) {
		case Element::CODE::SQUARE4:
			return 9;
			break;
		case Element::CODE::SQUARE8:
			return 23;
			break;
		case Element::CODE::TRIANGLE3:
			return 5;
			break;
		case Element::CODE::TRIANGLE6:
			return 22;
			break;
		case Element::CODE::TETRA4:
			return 10;
			break;
		case Element::CODE::TETRA10:
			return 24;
			break;
		case Element::CODE::PYRAMID5:
			return 14;
			break;
		case Element::CODE::PYRAMID13:
			return 27;
			break;
		case Element::CODE::PRISMA6:
			return 13;
			break;
		case Element::CODE::PRISMA15:
			return 26;
			break;
		case Element::CODE::HEXA8:
			return 12;
			break;
		case Element::CODE::HEXA20:
			return 25;
			break;
		default:
			return -1;
		}
	}
};
}



#endif /* SRC_OUTPUT_RESULT_VISUALIZATION_VTKWRITTER_H_ */
