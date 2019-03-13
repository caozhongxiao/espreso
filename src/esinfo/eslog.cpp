
#include "eslog.hpp"
#include "basis/logging/timelogger.h"
#include "basis/logging/progresslogger.h"
#include "basis/logging/solverlogger.h"
#include "basis/logging/oldtimelogger.h"
#include "basis/utilities/parser.h"
#include "basis/utilities/sysutils.h"
#include "basis/utilities/communication.h"
#include "config/reader/tokenizer.h"

#include "esinfo/mpiinfo.h"
#include "basis/logging/logger.h"

#include <unistd.h>
#include <ctime>
#include <cstring>
#include <string>
#include <fstream>

namespace espreso{
namespace eslog {

struct LoggerData {
	time_t initTime;
	std::string ecf;
	std::string outputRoot;
	std::string outputDirectory; // datetime
	std::string outputPath; // root/directory

	std::string name; // ecf without suffix
	std::string logFile; // root/directory/name.log
};

struct Logger: public espreso::Logger<OldTimeLogger, TimeLogger, ProgressLogger>, public LoggerData, public SolverLogger {};
Logger *logger = NULL;

const char* path()
{
	return logger->outputPath.c_str();
}

const char* name()
{
	return logger->name.c_str();
}

double time()
{
	return logger->time();
}

double duration()
{
	return logger->duration();
}

bool printtime()
{
	return logger->OldTimeLogger::verbosity > 1;
}

void create()
{
	logger = new Logger();
	logger->initTime = std::time(NULL);
	logger->ecf = "espreso.ecf";
	logger->outputRoot = "results";
}

void init(int *argc, char ***argv)
{
	logger->rank = info::mpi::rank;
	logger->size = info::mpi::size;

	if (info::mpi::rank == 0) {
		int option;
		while ((option = getopt(*argc, *argv, ":c:")) != -1) {
			if (option == 'c') {
				logger->ecf = optarg;
				if (std::ifstream(optarg).good()) {
					Tokenizer tok(logger->ecf);
					bool output = false, path = false, inoutput = false, read = true;
					while (read) {
						switch (tok.next()) {
						case Tokenizer::Token::STRING:
							if (path) {
								logger->outputRoot = tok.value();
								read = false;
							}
							if (!output) {
								output = StringCompare::caseInsensitiveEq(tok.value(), "OUTPUT");
							} else {
								if (inoutput) {
									path = StringCompare::caseInsensitiveEq(tok.value(), "PATH");
								}
							}
							break;
						case Tokenizer::Token::OBJECT_OPEN:
							inoutput = output;
							break;
						case Tokenizer::Token::END:
							read = false;
							break;
						default:
							break;
						}
					}
				}
			}
		}
	}

	auto synchronize = [] (std::string &str) {
		size_t ssize = str.size();
		MPI_Bcast(&ssize, sizeof(size_t), MPI_BYTE, 0, info::mpi::comm);
		char* dir = new char[ssize];
		if (info::mpi::rank == 0) {
			std::memcpy(dir, str.c_str(), str.size());
		}
		MPI_Bcast(dir, ssize, MPI_CHAR, 0, info::mpi::comm);
		str = std::string(dir, dir + ssize);
		delete[] dir;
	};

	// synchronize data accross MPI ranks
	MPI_Bcast(&logger->initTime, sizeof(time_t), MPI_BYTE, 0, info::mpi::comm);
	synchronize(logger->ecf);
	synchronize(logger->outputRoot);

	// compute output directory (path + datatime)
	struct tm *timeinfo;
	char buf[80];
	timeinfo = std::localtime(&logger->initTime);
	std::strftime(buf, 80, "%F-at-%Hh-%Mm-%Ss", timeinfo);
	logger->outputDirectory = std::string(buf);
	logger->outputPath = logger->outputRoot + "/" + logger->outputDirectory;
	size_t namebegin = logger->ecf.find_last_of("/") + 1;
	size_t nameend = logger->ecf.find_last_of(".");
	logger->name = logger->ecf.substr(namebegin, nameend - namebegin);
	logger->logFile = logger->outputPath + "/" + logger->name + ".log";

	if (info::mpi::rank) {
		MPI_Barrier(info::mpi::comm);
		logger->setLogFile(logger->logFile.c_str());
		return;
	} else {
		utils::createDirectory(logger->outputPath);
		logger->setLogFile(logger->logFile.c_str());
		MPI_Barrier(info::mpi::comm);
	}

	std::string symlink = logger->outputRoot + "/last";
	if (utils::exists(symlink)) {
		utils::remove(symlink);
	}
	utils::createSymlink(logger->outputDirectory, symlink);
	utils::copyFile(logger->ecf, symlink + "/" + logger->name + ".ecf");
}

void finish()
{
	logger->evaluate(*logger);
	if (logger) delete logger;
}

void start(const char* name, const char* section)
{
	logger->start(name, section);
}

void checkpoint(const char* name)
{
	logger->checkpoint(name);
}

void end(const char* name)
{
	logger->end(name);
}

void ln()
{
	logger->ln();
}

void nextStep(int step)
{
	logger->nextLoadStep(step);
}

void startln(const char* name, const char* section)
{
	start(name, section);
	ln();
}

void checkpointln(const char* name)
{
	checkpoint(name);
	ln();
}

void endln(const char* name)
{
	end(name);
	logger->ln();
}

void param(const char* name, const int &value)
{
	logger->param(name, value);
}

void param(const char* name, const long &value)
{
	logger->param(name, value);
}

void param(const char* name, const long unsigned int &value)
{
	logger->param(name, value);
}

void param(const char* name, const double &value)
{
	logger->param(name, value);
}

void addsolverparam(const char* name, const char* shortcut, const char* format, int &value)
{
	logger->addparam(name, shortcut, format, value);
}

void addsolverparam(const char* name, const char* shortcut, const char* format, long &value)
{
	logger->addparam(name, shortcut, format, value);
}

void addsolverparam(const char* name, const char* shortcut, const char* format, long unsigned int &value)
{
	logger->addparam(name, shortcut, format, value);
}

void addsolverparam(const char* name, const char* shortcut, const char* format, double &value)
{
	logger->addparam(name, shortcut, format, value);
}

void addsolverparam(const char* name, const char* shortcut, const char* format, const char* &value)
{
	logger->addparam(name, shortcut, format, value);
}

void addsolverparam(const char* name, const char* shortcut, bool &value)
{
	logger->addparam(name, shortcut, value);
}

void printsolverheader()
{
	logger->printheader(progress());
}

void printsolver()
{
	logger->print(progress());
}

void param(const char* name, const char* value)
{
	logger->param(name, value);
}

void info(const char* msg)
{
	logger->info(msg);
}

void solver(const char* msg)
{
	logger->info(msg);
}

void linearsolver(const char* msg)
{
	if (progress().verbosity > 2) {
		logger->info(msg);
	}
}

void duration(const char* msg)
{
	if (printtime()) {
		logger->info(msg);
	}
}

void warning(const char* msg)
{
	logger->warning(msg);
}

void debug(const char* msg)
{
	logger->debug(msg);
}

void error(const char* msg)
{
	logger->error(msg);
}

void globalerror(const char* msg)
{
	if (logger->rank == 0) {
		logger->error(msg);
	}
}

ProgressLogger& progress()
{
	return *logger;
}

}
}



