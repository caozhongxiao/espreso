
#include "sysutils.h"

#include "esinfo/mpiinfo.h"
#include "esinfo/timeinfo.h"
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

	createDirectory(prefix.str());
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

std::string debugDirectory()
{
	std::stringstream path;
	path << eslog::path() << "/DEBUG";
	path << "/step" << time::step;
	path << "/substep" << time::substep;
	path << "/iteration" << time::iteration;
	path << "/" << info::mpi::rank;
	return path.str();
}

std::string prepareFile(const std::string &directory, const std::string &name, int domain)
{
	createDirectory(directory);
	if (domain != -1) {
		return std::string(directory + "/" + name + std::to_string(domain) + ".txt");
	} else {
		return std::string(directory + "/" + name + ".txt");
	}
}

void printStack()
{
	pid_t pid = getpid();
	std::string pstack = "pstack " + std::to_string(pid);
	system(pstack.c_str());
}

}
}





