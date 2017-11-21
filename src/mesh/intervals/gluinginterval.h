
#ifndef SRC_MESH_INTERVALS_GLUINGINTERVAL_H_
#define SRC_MESH_INTERVALS_GLUINGINTERVAL_H_

#include "domaininterval.h"

namespace espreso {

struct GluingInterval: public DomainInterval {
	eslocal dindex, ndomains;
	eslocal LMOffset;

	GluingInterval(const DomainInterval &interval, eslocal dindex, eslocal ndomains)
	: DomainInterval(interval), dindex(dindex), ndomains(ndomains), LMOffset(-1) {}
};

}


#endif /* SRC_MESH_INTERVALS_GLUINGINTERVAL_H_ */
