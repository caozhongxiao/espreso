
#include "sysutils.h"

#include "esinfo/eslog.hpp"

#include <execinfo.h>
#include <cstdlib>
#include <unistd.h>
#include <sstream>
#include <fstream>
#include <algorithm>

namespace espreso {
namespace utils {

std::string createDirectory(const std::vector<std::string> &path)
{
	std::stringstream prefix;
	std::for_each(path.begin(), path.end(), [&] (const std::string &dir) { prefix << dir << "/"; });

	if (std::system(("mkdir -p " + prefix.str()).c_str())) {
		eslog::error("Cannot create directory '%s'\n", prefix.str().c_str());
	}
	return prefix.str();
}

void createDirectory(const std::string &path)
{
	if (std::system(("mkdir -p " + path).c_str())) {
		eslog::error("Cannot create directory '%s'\n", path.c_str());
	}
}

bool exists(const std::string &path)
{
	return std::ifstream(path.c_str()).good();
}

void remove(const std::string &path)
{
	if (std::remove(path.c_str())) {
		eslog::error("Cannot remove '%s'\n", path.c_str());
	}
}

void createSymlink(const std::string &path, const std::string &link)
{
	if (symlink(path.c_str(), link.c_str())) {
		eslog::error("Cannot create link '%s' to '%s'\n", link.c_str(), path.c_str());
	}
}

void copyFile(const std::string &source, const std::string &destination)
{
	std::ifstream src(source.c_str(), std::ios::binary);
	if (!src.good()) {
		eslog::error("Cannot read file '%s'\n", source.c_str());
	}
	std::ofstream dst(destination, std::ios::binary);
	if (!dst.good()) {
		eslog::error("Cannot create file '%s'\n", destination.c_str());
	}
	dst << src.rdbuf();
}

void currentStack(std::vector<std::string> &stack)
{
	stack.reserve(30);
	std::vector<void*> stackdata(30);
	size_t size = backtrace(stackdata.data(), 30);
	char** functions = backtrace_symbols(stackdata.data(), size);

	std::stringstream command;
	command << "addr2line -sipfC -e $(which espreso)";
	for (size_t i = 0; i < size; i++) {
		std::string function(functions[i]);
		size_t begin = function.find_last_of('[') + 1;
		size_t end = function.find_last_of(']');
		command << " " << function.substr(begin, end - begin);
	}
	free(functions);

	FILE *in;
	char buff[512];
	if(!(in = popen(command.str().c_str(), "r"))){
		return; // broken
	}

	std::string message;
	while(fgets(buff, sizeof(buff), in) != NULL){
		stack.push_back(buff);
	}
	pclose(in);
}

}
}





