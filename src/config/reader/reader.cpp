
#include "reader.h"
#include "tokenizer.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "config/ecf/output.h"
#include "config/configuration.h"
#include "config/holders/valueholder.h"
#include "basis/logging/verbosity.h"
#include "basis/utilities/utils.h"
#include "basis/utilities/parser.h"

#include "mpi.h"

#include <getopt.h>
#include <stack>
#include <functional>
#include <unistd.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstring>

using namespace espreso;

template<typename T>
std::ostream& operator<< (std::ostream& os, const std::vector<T> &v)
{
	if (v.size()) {
		os << v[0];
	}
	for(size_t i = 1; i < v.size(); ++i) {
		os << "::" << v[i];
	}
	return os;
}

std::string ECFReader::configurationFile = "espreso.ecf";
std::vector<VerboseArg*> ECFReader::verbosity;

static struct option long_options[] = {
		{"config",  required_argument, 0, 'c'},
		{"default",  optional_argument, 0, 'd'},
		{"help",  no_argument, 0, 'h'},
		{0, 0, 0, 0}
};

ECFReader::ECFReader(int argc, char **argv)
: ecf("espreso.ecf"), outputPath("results"), meshDuplication(1)
{
	if (info::mpi::rank == 0) {
		int option;
		while ((option = getopt(argc, argv, ":c:")) != -1) {
			if (option == 'c') {
				ecf = optarg;
				if (std::ifstream(optarg).good()) {
					Tokenizer tok(ecf);
					bool output = false, path = false, inoutput = false, read = true;
					bool decomposition = false, duplication = false, indecomposition = false;
					while (read) {
						switch (tok.next()) {
						case Tokenizer::Token::STRING:
							if (path) {
								outputPath = tok.value();
								path = false;
							}
							if (duplication) {
								ECFValueHolder<int> dup(meshDuplication);
								if (!dup.setValue(tok.value())) {
									eslog::error("Invalid duplication value (not int).\n");
								}
								duplication = false;
							}
							if (!output) {
								output = StringCompare::caseInsensitiveEq(tok.value(), "OUTPUT");
							} else {
								if (inoutput) {
									path = StringCompare::caseInsensitiveEq(tok.value(), "PATH");
								}
							}
							if (!decomposition) {
								decomposition = StringCompare::caseInsensitiveEq(tok.value(), "DECOMPOSITION");
							} else {
								if (indecomposition) {
									duplication = StringCompare::caseInsensitiveEq(tok.value(), "MESH_DUPLICATION");
								}
							}
							break;
						case Tokenizer::Token::OBJECT_OPEN:
							inoutput = output;
							indecomposition = decomposition;
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
	MPI_Bcast(&meshDuplication, sizeof(int), MPI_BYTE, 0, info::mpi::comm);
	synchronize(ecf);
	synchronize(outputPath);
}

ECFReader::ECFReader(const std::string &file)
: ecf(file), outputPath("results"), meshDuplication(1)
{

}

ECFRedParameters ECFReader::_read(
		ECFObject &configuration,
		int* argc,
		char ***argv,
		const std::map<size_t, std::string> &defaultArgs,
		const std::map<std::string, std::string> &variables)
{
	int option_index, option;
	std::string options("c:dhvm");

	std::vector<struct option> opts;
	std::vector<std::pair<std::string, ECFParameter*> > parameters;
	std::vector<std::string> nameless;

	std::function<void(const ECFObject &configuration, std::vector<std::string> path)>
	recurse = [&] (const ECFObject &configuration, std::vector<std::string> path) {
		for (size_t i = 0; i < configuration.parameters.size(); i++) {
			if (configuration.parameters[i]->isValue()) {
				std::string prefix;
				std::for_each(path.begin(), path.end(), [&] (const std::string &p) { prefix += p + "::"; });
				parameters.push_back(std::make_pair(prefix + Parser::uppercase(configuration.parameters[i]->name), configuration.parameters[i]));
			}
			if (configuration.parameters[i]->isObject()) {
				path.push_back(Parser::uppercase(configuration.parameters[i]->name));
				recurse(dynamic_cast<const ECFObject&>(*configuration.parameters[i]), path);
				path.pop_back();
			}
		}
	};
	recurse(configuration, {});

	opts.reserve(parameters.size() + 3);
	for (size_t i = 0; i < parameters.size(); i++) {
		opts.push_back({ parameters[i].first.c_str(), required_argument, 0, 'p' });
	}

	option_index = 0;
	while (long_options[option_index].name != 0) {
		opts.push_back(long_options[option_index++]);
	}

	// read the rest parameters
	size_t helpVerboseLevel = 0;

	std::string confFile = "espreso.ecf";
	std::string sopts = "c:d::h";

	for (size_t i = 0; i < verbosity.size(); i++) {
		sopts.append({ verbosity[i]->argflag });
	}

	optind = 0;
	while ((option = getopt_long(*argc, *argv, sopts.c_str(), opts.data(), &option_index)) != -1) {
		switch (option) {
		case 'p':
			// parameters will be read after configuration file
			break;
		case 'h':
			helpVerboseLevel++;
			break;
		case 'd': {
			std::ofstream os("espreso.ecf.default");
			store(configuration, os, false, true);
			os.close();
			exit(EXIT_SUCCESS);
		} break;
		case 'c':
			confFile = optarg;
			break;
		case '?':
			exit(EXIT_FAILURE);
			break;
		}
	}

	if (helpVerboseLevel) {
		std::cout << "\nusage: espreso [options] [ARGS]\n\n";

		std::cout << " [options] are the following:\n";
		std::cout << "\t -h, --help            print this message\n";
		std::cout << "\t -d, --default         generate default configuration file\n";
		std::cout << "\t -c, --config=[path]   set path to configuration file\n\n";

		std::cout << " [ARGS] are unnamed argument that can be referenced in a configuration\n";
		std::cout << "        file by [ARGX], where X is a number of an argument counted from 0.\n";
		std::cout << "        e.g. DOMAINS [ARG0]; set number of domains to the first argument.\n\n";

		std::cout << "The solver is controlled by '*.ecf' scripts. Some examples can be found\n";
		std::cout << "in 'benchmarks' directory or in 'tests/examples'. In default 'espreso.ecf'\n";
		std::cout << "from ESPRESO root directory is loaded. A different configuration file can\n";
		std::cout << "be set by -c [path] option.\n\n";
		std::cout << "A file is composed from parameters and objects with the following pattern:\n\n";
		std::cout << "  OBJECT {\n";
		std::cout << "    PARAMETER VALUE;\n";
		std::cout << "  }\n\n";
		std::cout << "Run espreso --default to generates the configuration file with default\n";
		std::cout << "parameters.\n\n";
		exit(0);
	}

	// read nameless parameters
	while (optind < *argc) {
		nameless.push_back(std::string((*argv)[optind++]));
	}
	configurationFile = confFile;

	ECFRedParameters redParameters = _read(configuration, confFile, nameless, defaultArgs, variables);

	optind = 0;
	while ((option = getopt_long(*argc, *argv, sopts.c_str(), opts.data(), &option_index)) != -1) {
		switch (option) {
		case 'p':
			if (!parameters[option_index].second->setValue(optarg)) {
				eslog::globalerror("Parameter '%s' has wrong value '%s'\n", parameters[option_index].first.c_str(), optarg);
			}
			break;
		case 'h':
		case 'd':
		case 'c': break;
		default:
			for (size_t i = 0; i < verbosity.size(); i++) {
				if (verbosity[i]->argflag == option) {
					++verbosity[i]->verbosity;
					break;
				}
			}
		}
	}

	return redParameters;
}

ECFRedParameters ECFReader::_read(
		ECFObject &configuration,
		const std::string &file,
		const std::vector<std::string> &args,
		const std::map<size_t, std::string> &defaultArgs,
		const std::map<std::string, std::string> &variables)
{
	ECFRedParameters redParameters;
	if (info::mpi::rank == 0) {
		std::ifstream ecffile(file);
		redParameters.hadValidECF = ecffile.good();
	}
	int valid = redParameters.hadValidECF;
	MPI_Bcast(&valid, 1, MPI_INT, 0, info::mpi::comm);
	redParameters.hadValidECF = valid;

	if (!redParameters.hadValidECF) {
		return redParameters;
	}

	std::vector<std::string> prefix;
	std::vector<std::string> values;
	std::stack<ECFObject*> confStack;
	std::stack<Tokenizer*> tokenStack;

	bool correctlyLoaded = true, unfinishedLine = false;
	std::map<size_t, std::vector<std::string> > arguments;

	std::string link;
	ECFParameter *parameter;
	confStack.push(&configuration);
	tokenStack.push(new CollectiveTokenizer(file));
	while (tokenStack.size()) {
		switch (tokenStack.top()->next()) {
		case Tokenizer::Token::END:
			delete tokenStack.top();
			tokenStack.pop();
			break;
		case Tokenizer::Token::STRING:
			values.push_back(tokenStack.top()->value());
			if (unfinishedLine && confStack.top()->getParameter(values.back()) != NULL) {
				eslog::globalerror("PARSE ERROR: Expected ';' before '%s'\n%s", values.back().c_str(), tokenStack.top()->lastLines(2).c_str());
			}
			break;
		case Tokenizer::Token::LINK:
		{
			std::string value = tokenStack.top()->value();
			link = "[" + value + "]";
			for (auto it = variables.begin(); it != variables.end(); ++it) {
				if (StringCompare::caseInsensitiveEq(it->first, value)) {
					value = it->second;
					break;
				}
			}
			if (value.size() > 4 && StringCompare::caseInsensitiveSuffix(value, ".ecf")) {
				tokenStack.push(new CollectiveTokenizer(value));
				break;
			}
			if (value.size() > 2 && StringCompare::caseInsensitivePreffix("ARG", value)) {
				std::stringstream ss(std::string(value.begin() + 3, value.end()));
				size_t index;
				ss >> index;
				if (!ss.fail() && ss.eof() && index < args.size()) {
					values.push_back(args[index]);
					std::string parameter;
					std::for_each(prefix.begin(), prefix.end(), [&] (const std::string &s) { parameter += s + "::"; });
					arguments[index].push_back(parameter + values.front());
				} else {
					if (index < args.size()) {
						eslog::globalerror("Invalid argument '%s'\n", value.c_str());
					} else {
						auto ait = defaultArgs.find(index);
						if (ait != defaultArgs.end()) {
							values.push_back(ait->second);
						} else {
							correctlyLoaded = false;
							if (values.size()) {
								std::string parameter;
								std::for_each(prefix.begin(), prefix.end(), [&] (const std::string &s) { parameter += s + "::"; });
								arguments[index].push_back(parameter + values.front());
							} else {
								eslog::globalerror("Parameter cannot be the [ARG].\n", tokenStack.top()->lastLines(2).c_str());
							}
						}
					}
				}
				break;
			}
			if (value.size() > 4 && StringCompare::caseInsensitiveSuffix(value, ".csv")) {
				CollectiveTokenizer csv(value);
				bool run = true;
				while(run) {
					switch (csv.next()) {
					case Tokenizer::Token::END:
						run = false;
						break;
					case Tokenizer::Token::STRING:
						values.push_back(csv.value() + (values.size() % 2 == 0 ? "," : ";"));
						break;
					case Tokenizer::Token::LINE_END:
					case Tokenizer::Token::DELIMITER:
					case Tokenizer::Token::EXPRESSION_END:
						break;
					default:
						eslog::globalerror("Error while reading file '%s'\n", value.c_str());
					}
				}
				break;
			}
			if (StringCompare::caseInsensitiveEq(values.back(), "TABULAR")) {
				link = "TABULAR [" + value + "]";
				values.push_back("[" + value + "]");
			} else {
				values.push_back(value);
			}
			break;
		}
		case Tokenizer::Token::OBJECT_OPEN:
			if (values.size() == 0) {
				eslog::globalerror("PARSE ERROR: Opening of an unnamed region is not allowed.\n%s", tokenStack.top()->lastLines(2).c_str());
			}
			if (values.size() > 1) {
				eslog::globalerror("PARSE ERROR: Multiple names for a region are not allowed.\n%s", tokenStack.top()->lastLines(2).c_str());
			}
			prefix.push_back(values[0]);
			parameter = confStack.top()->getParameter(values[0]);
			if (parameter == NULL) {
				std::string prefixstr;
				std::for_each(prefix.begin(), prefix.end(), [&] (const std::string &s) { prefixstr += s + "::"; });
				eslog::globalerror("PARSE ERROR: Unexpected parameter '%s'\n%s", prefixstr.c_str(), tokenStack.top()->lastLines(2).c_str());
			}
			if (!parameter->isObject()) {
				std::string prefixstr;
				std::for_each(prefix.begin(), prefix.end(), [&] (const std::string &s) { prefixstr += s + "::"; });
				eslog::globalerror("PARSE ERROR: Expected parameter instead of object '%s'.\n%s", prefixstr.c_str(), tokenStack.top()->lastLines(2).c_str());
			}
			redParameters.parameters.push_back(parameter);
			confStack.push(dynamic_cast<ECFObject*>(parameter));
			values.clear();
			break;
		case Tokenizer::Token::OBJECT_CLOSE:
			if (!confStack.size() || !prefix.size()) {
				eslog::globalerror("PARSE ERROR: Unexpected region end.\n%s", tokenStack.top()->lastLines(2).c_str());
			}
			prefix.pop_back();
			confStack.pop();
			break;
		case Tokenizer::Token::ASSIGN:
			break;
		case Tokenizer::Token::DELIMITER:
			values.push_back(",");
			break;
		case Tokenizer::Token::EXPRESSION_END:
		{
			unfinishedLine = false;
			if (!correctlyLoaded) {
				values.clear();
				break;
			}
			if (values.size() == 0) {
				break;
			}
			if (values.size() == 1) {
				// allow to read an empty value
				values.push_back("");
			}
			if (values.size() < 2) {
				eslog::globalerror(
						"PARSE ERROR: Incorrect assignment format on line %ld. Use 'PARAMETER' 'VALUE';\n%s",
						tokenStack.top()->line(), tokenStack.top()->lastLines(2).c_str());
			}
			std::stringstream ss;
			ss << values[1];
			for (size_t i = 2; i < values.size(); i++) {
				ss << " " << values[i];
			}
			parameter = confStack.top()->getParameter(values[0]);
			if (parameter == NULL) {
				prefix.push_back(values[0]);
				std::stringstream ssp; ssp << prefix;
				eslog::globalerror("PARSE ERROR: Unexpected parameter '%s'\n%s", ssp.str().c_str(), tokenStack.top()->lastLines(2).c_str());
			}
			if (parameter->isObject()) {
				std::stringstream ssp; ssp << prefix;
				eslog::globalerror("Invalid ECF configuration. Parameter '%s::%s' is an object. Expected '{' instead of '%s'\n",
						ssp.str().c_str(), parameter->name.c_str(), tokenStack.top()->lastLines(2).c_str(), ss.str().c_str());
			}
			if (!parameter->setValue(ss.str())) {
				eslog::globalerror("PARSE ERROR: Parameter '%s' has wrong value '%s'\n", values[0].c_str(), ss.str().c_str());
			}
			redParameters.parameters.push_back(parameter);
			if (link.size()) {
				redParameters.defaulted[redParameters.parameters.back()] = link;
				link = "";
			}
			values.clear();
			break;
		}
		case Tokenizer::Token::LINE_END:
			if (values.size() > 1) {
				unfinishedLine = true;
			}
			break;
		default:
			eslog::globalerror("PARSE ERROR: Unknown token in configuration file.");
		}
	}
	if (confStack.size() != 1) {
		eslog::globalerror("PARSE ERROR: Unexpected EOF before close all regions.");
	}

	if (!correctlyLoaded) {
		std::string error = "Configuration file is not correctly loaded.\nUse espreso ";
		size_t i = 0;
		for (auto it = arguments.begin(); it != arguments.end(); ++it) {
			error += "[ARG" + std::to_string(i++) + "] ";
		}
		error += "\nWhere ARGs are the following:\n";
		i = 0;
		for (auto it = arguments.begin(); it != arguments.end(); ++it) {
			error += "ARG" + std::to_string(i++) + " = { ";
			for (size_t j = 0; j < it->second.size(); j++) {
				error += it->second[j];
				if (j != it->second.size() - 1) {
					error += ", ";
				}
			}
			error += " }\n";
		}
		eslog::globalerror(error.c_str());
	}
	return redParameters;
}

static void printECF(const ECFObject &configuration, std::ostream &os, size_t indent, bool hasDataType, bool onlyAllowed, bool printPatterns, bool pattern, const ECFRedParameters &parameters)
{
	auto printindent = [&] (size_t indent) {
		for (size_t j = 0; j < indent; j++) {
			os << " ";
		}
	};

	auto wasRed = [&] (const ECFParameter* parameter) {
		if (parameters.parameters.size()) {
			return std::find(parameters.parameters.begin(), parameters.parameters.end(), parameter) != parameters.parameters.end();
		}
		return true;
	};

	size_t maxSize = 0;
	for (size_t i = 0; i < configuration.parameters.size(); i++) {
		if (!wasRed(configuration.parameters[i])) {
			continue;
		}
		if ((configuration.parameters[i]->metadata.isallowed() || !onlyAllowed) && configuration.parameters[i]->isValue()) {
			std::string value = configuration.parameters[i]->getValue();
			if (parameters.defaulted.find(configuration.parameters[i]) != parameters.defaulted.end()) {
				value = parameters.defaulted.find(configuration.parameters[i])->second;
			}
			if (maxSize < configuration.parameters[i]->name.size() + value.size()) {
				maxSize = configuration.parameters[i]->name.size() + value.size();
			}
		}
	}

	bool printSpace = false;
	bool firstParameter = true;
	bool firstValue = true;

	auto printparameter = [&] (const ECFParameter *parameter) {
		if (onlyAllowed && !parameter->metadata.isallowed()) {
			return;
		}
		if (!wasRed(parameter)) {
			// Empty spaces are never red. Hence, they cannot be skipped
			if (parameter->isObject() || parameter->isValue()) {
				return;
			}
		}

		if (parameter->isValue()) {
			if (printSpace) {
				if (!firstValue) {
					os << "\n";
				}
				printSpace = false;
			}
			std::string value = parameter->getValue();
			if (parameters.defaulted.find(parameter) != parameters.defaulted.end()) {
				value = parameters.defaulted.find(parameter)->second;
			}
			size_t space = maxSize ? maxSize - parameter->name.size() - value.size() : 3;
			if (printPatterns && parameter->metadata.datatype.front() == ECFDataType::OPTION) {
				printindent(indent);
				os << "#[";
				for (size_t i = 0; i < parameter->metadata.options.size(); i++) {
					if (i) {
						os << ",";
					}
					os << Parser::uppercase(parameter->metadata.options[i].name);
				}
				os << "]\n";
			}
			printindent(indent);
			if (hasDataType) {
				os << parameter->name;
			} else {
				os << Parser::uppercase(parameter->name);
			}
			printindent(space + 3);
			if (
					parameter->metadata.datatype.front() == ECFDataType::STRING ||
					parameter->metadata.datatype.front() == ECFDataType::BOUNDARY_REGION ||
					parameter->metadata.datatype.front() == ECFDataType::MATERIAL) {
				os << value << ";\n";
			} else {
				os << Parser::uppercase(value) << ";\n";
			}
			firstValue = false;
		} else if (parameter->isObject()) {
			if (!firstParameter) {
				os << "\n";
			}
			printSpace = false;
			printindent(indent);
			if (hasDataType) {
				os << parameter->name << " {\n";
			} else {
				os << Parser::uppercase(parameter->name) << " {\n";
			}
			printECF(*dynamic_cast<const ECFObject*>(parameter), os, indent + 2, parameter->metadata.datatype.size(), onlyAllowed, printPatterns, pattern, parameters);
			printindent(indent);
			os << "}\n";
			firstValue = false;
		} else {
			printSpace = true;
			// Separators etc..
			return;
		}
		firstParameter = false;
	};

	bool spaceAfterObject = false;
	for (size_t i = 0; i < configuration.parameters.size(); i++) {
		if (!configuration.parameters[i]->isObject() && !configuration.parameters[i]->isValue()) {
			spaceAfterObject = false;
		}
		if (spaceAfterObject && (wasRed(configuration.parameters[i]))) {
			spaceAfterObject = false;
			if (configuration.parameters[i]->isValue()) {
				os << "\n";
			}
		}
		if (configuration.parameters[i]->isObject() && wasRed(configuration.parameters[i])) {
			spaceAfterObject = true;
		}
		printparameter(configuration.parameters[i]);
	}
	if (printPatterns && configuration.getPattern()) {
		if (pattern == false) {
			pattern = true;
			printindent(indent);
			os << "/*\n";
			printparameter(configuration.getPattern());
			printindent(indent);
			os << "*/\n";
		} else {
			printparameter(configuration.getPattern());
		}
	}
};

void ECFReader::store(const ECFObject &configuration, std::ostream &os, bool onlyAllowed, bool printPatterns, const ECFRedParameters &parameters)
{
	os << "# ESPRESO Configuration File\n\n";
	printECF(configuration, os, 0, false, onlyAllowed, printPatterns, false, parameters);
}
