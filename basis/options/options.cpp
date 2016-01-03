
#include "options.h"

static struct option long_options[] = {
		{"input",  required_argument, 0, 'i'},
		{"path",  required_argument, 0, 'p'},
		{"help",  no_argument, 0, 'h'},
		{0, 0, 0, 0}
};

Options::Options(int* argc, char*** argv): verbosity(VERBOSE)
{
	auto printOption = [] (const std::string &opt, const std::string &desc) {
		std::cout << "\t" << opt;
		if (opt.length() < 16) {
			std::cout << "\t\t" << desc << "\n";
		} else {
			std::cout << "\t" << desc << "\n";
		}
	};

	int option_index, option;
	while (true) {
		option = getopt_long(*argc, *argv, "i:p:vh", long_options, &option_index);
		if (option == -1) {
			break;
		}

		switch (option) {
		case 'v':
			verbosity++;
			break;
		case 'i':
			input = std::string(optarg);
			input.erase(0, input.find_first_not_of('='));
			break;
		case 'p':
			path = std::string(optarg);
			path.erase(0, path.find_first_not_of('='));
			break;
		case 'h':
			std::cout << "Usage: espreso [OPTIONS] [PARAMETERS]\n";

			std::cout << "\nOPTIONS:\n";
			printOption("-h, --help", "show this message");
			printOption("-i, --input=INPUT", "input format: [generator, matsol, workbench, esdata]");
			printOption("-p, --path=PATH", "path to an example");
			printOption("-v,vv,vvv", "verbose level");

			std::cout << "\nPARAMETERS:\n";
			std::cout << "\tlist of nameless parameters for particular example\n";
			break;
		case '?':
			break;
		}
	}

	while (optind < *argc) {
		nameless.push_back(std::string((*argv)[optind++]));
	}
}

std::ostream& operator<<(std::ostream& os, const Options &options)
{
       os << "input: '" << options.input << "'\n";
       os << "path: '" << options.path << "'\n";
       os << "verbosity level: " << options.verbosity << "\n";
       return os;
}




