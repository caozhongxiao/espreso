
#include "logging.h"

namespace espreso {

Test::~Test()
{
	if (!error) {
		return;
	}

	ESINFO(ERROR) << "ESPRESO INTERNAL TEST FAILED: "<< os.str();
}


static void printStack()
{
	std::vector<void*> stack(30);
	size_t size = backtrace(stack.data(), 30);
	char** functions = backtrace_symbols(stack.data(), size);

	std::stringstream command;
	command << "addr2line -sipfC -e " << config::env::executable;
	for (size_t i = 0; i < size; i++) {
		std::string function(functions[i]);
		size_t begin = function.find_last_of('[') + 1;
		size_t end = function.find_last_of(']');
		command << " " << function.substr(begin, end - begin);
	}
	free(functions);
	if (system(command.str().c_str())) { // convert addresses to file lines
		ESINFO(ALWAYS) << "Broken address to file lines command";
	}
}


Info::~Info()
{
	if (_plain) {
		if (config::env::MPIrank == 0) {
			std::cout << os.str();
		}
		return;
	}
	if (event == ERROR || (event == GLOBAL_ERROR && config::env::MPIrank == 0)) {
		fprintf(stderr, "\x1b[31m%s\x1b[0m\n", os.str().c_str());
		if (event == ERROR) {
			fprintf(stderr, "ESPRESO EXITED WITH AN ERROR ON PROCESS %d.\n\n\n", config::env::MPIrank);

			if (config::env::executable.size()) {
				printStack();
			}
		}

		fflush(stderr);
		exit(EXIT_FAILURE);
	}

	os << std::endl;

	if (config::env::MPIrank != 0) {
		return; // only first process print results
	}

	switch (color) {
	case TextColor::WHITE:
		fprintf(stdout, "%s", os.str().c_str());
		break;
	case TextColor::RED:
		fprintf(stdout, "\x1b[31m%s\x1b[0m", os.str().c_str());
		break;
	case TextColor::GREEN:
		fprintf(stdout, "\x1b[32m%s\x1b[0m", os.str().c_str());
		break;
	case TextColor::YELLOW:
		fprintf(stdout, "\x1b[33m%s\x1b[0m", os.str().c_str());
		break;
	case TextColor::BLUE:
		fprintf(stdout, "\x1b[34m%s\x1b[0m", os.str().c_str());
		break;
	}

	fflush(stdout);
}


Measure::~Measure()
{
//	switch (event) {
//	case CHECKPOINT3:
//		checkpoints.push_back(Checkpoint(os.str(), time(), 3));
//		break;
//	case CHECKPOINT2:
//		checkpoints.push_back(Checkpoint(os.str(), time(), 2));
//		break;
//	case CHECKPOINT1:
//		checkpoints.push_back(Checkpoint(os.str(), time(), 1));
//		break;
//	case SUMMARY:
//		checkpoints.push_back(Checkpoint(os.str(), time(), 0));
//		evaluateCheckpoints();
//		return;
//	}

	os << std::endl;

	if (config::env::MPIrank != 0) {
		return; // only first process print results
	}

	fprintf(stdout, "%s", os.str().c_str());
	fflush(stdout);
}

double Measure::processMemory()
{
	std::ifstream file("/proc/self/status");
	eslocal result = -1;
	std::string line, label("VmRSS:");

	while (getline(file, line)) {
		if (line.find(label.c_str(), 0, label.size()) == 0) {
			file.close();
			std::stringstream(line) >> label >> result;
			return result / 1024.0;
		}
	}
	file.close();
	return 0;
}

double Measure::usedRAM()
{
	struct sysinfo memInfo;
	sysinfo(&memInfo);
	return ((memInfo.totalram - memInfo.freeram) * memInfo.mem_unit) / 1024.0 / 1024.0;
}

double Measure::availableRAM()
{
	struct sysinfo memInfo;
	sysinfo(&memInfo);
	return (memInfo.totalram * memInfo.mem_unit) / 1024.0 / 1024.0;
}

}



