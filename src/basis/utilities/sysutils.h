
#ifndef SRC_BASIS_UTILITIES_SYSUTILS_H_
#define SRC_BASIS_UTILITIES_SYSUTILS_H_

#include <string>
#include <vector>

namespace espreso {
namespace utils {

std::string createDirectory(const std::vector<std::string> &path);
void createDirectory(const std::string &path);
void createSymlink(const std::string &path, const std::string &link);
void copyFile(const std::string &source, const std::string &destination);
bool exists(const std::string &path);
void remove(const std::string &path);

void printStack();

}
}


#endif /* SRC_BASIS_UTILITIES_SYSUTILS_H_ */
