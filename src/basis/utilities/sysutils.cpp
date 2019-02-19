
#include <basis/utilities/sysutils.h>
#include "basis/logging/logging.h"

#include <fstream>
#include <algorithm>

namespace espreso {
namespace utils {

std::string createDirectory(const std::vector<std::string> &path)
{
	std::stringstream prefix;
	std::for_each(path.begin(), path.end(), [&] (const std::string &dir) { prefix << dir << "/"; });

	if (system(("mkdir -p " + prefix.str()).c_str())) {
		ESINFO(ERROR) << "Cannot create requested directory";
	}
	return prefix.str();
}

}
}





