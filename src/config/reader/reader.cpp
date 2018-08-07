
#include "reader.h"

#include <getopt.h>
#include <stack>
#include <functional>
#include <regex>
#include <unistd.h>
#include <fstream>

#include "mpi.h"

#include "tokenizer.h"
#include "../ecf/output.h"
#include "../ecf/environment.h"
#include "../../basis/logging/logging.h"
#include "../../basis/utilities/parser.h"

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

static struct option long_options[] = {
		{"config",  required_argument, 0, 'c'},
		{"default",  optional_argument, 0, 'd'},
		{"help",  no_argument, 0, 'h'},
		{0, 0, 0, 0}
};

ECFRedParameters ECFReader::_read(
		ECFObject &configuration,
		int* argc,
		char ***argv,
		const std::map<size_t, std::string> &defaultArgs,
		const std::map<std::string, std::string> &variables)
{
	environment->executable = *argv[0];
	int option_index, option;
	std::string options("c:dhvtm");

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
	if (StringCompare::caseSensitiveSuffix(std::string(*argv[0]), "espreso")) {
		confFile = "espreso.ecf";
	}
	if (StringCompare::caseSensitiveSuffix(std::string(*argv[0]), "decomposer")) {
		confFile = "decomposer.ecf";
	}

	optind = 0;
	while ((option = getopt_long(*argc, *argv, "c:d::hvtm", opts.data(), &option_index)) != -1) {
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

	size_t start = confFile.find_last_of("/") + 1;
	size_t end   = confFile.find_last_of(".");
	Logging::name = confFile.substr(start, end - start);
	configurationFile = confFile;

	ECFRedParameters redParameters = _read(configuration, confFile, nameless, defaultArgs, variables);

	optind = 0;
	while ((option = getopt_long(*argc, *argv, "c:dhvtm", opts.data(), &option_index)) != -1) {
		switch (option) {
		case 'p':
			if (!parameters[option_index].second->setValue(optarg)) {
				ESINFO(GLOBAL_ERROR) << "Parameter '" << parameters[option_index].first << "' has wrong value '" << optarg << "'";
			}
			break;
		case 'v':
			environment->verbose_level++;
			break;
		case 't':
			environment->testing_level++;
			break;
		case 'm':
			environment->measure_level++;
			break;
		}
	}

	return redParameters;
}

void ECFReader::copyInputData()
{
	if (environment->MPIrank) {
		MPI_Barrier(environment->MPICommunicator);
		Logging::log.open(Logging::outputRoot() + "/" + Logging::name + ".log", std::ofstream::app);
		return;
	} else {
		if (environment->remove_old_results && Logging::path.compare(".")) {
			system(("rm -fr " + Logging::path).c_str());
		}
		if (system(("mkdir -p " + Logging::outputRoot()).c_str())) {
			ESINFO(ERROR) << "Cannot create output directory\n";
		}
		MPI_Barrier(environment->MPICommunicator);
	}

	Logging::log.open(Logging::outputRoot() + "/" + Logging::name + ".log", std::ofstream::app);

	int error = remove(std::string(Logging::path + "/" + "last").c_str());
	error = symlink(("../" + Logging::outputRoot()).c_str(), std::string(Logging::path + "/" + "last").c_str());
	if (error) {
		ESINFO(ALWAYS_ON_ROOT) << Info::TextColor::YELLOW << "Something wrong happens with creating link to last output directory.";
	}

	std::ifstream src(configurationFile.c_str(), std::ios::binary);
	std::ofstream dst((Logging::outputRoot() + "/" + configurationFile.substr(configurationFile.find_last_of("/") + 1)).c_str(), std::ios::binary);

	dst << src.rdbuf();
}

ECFRedParameters ECFReader::_read(
		ECFObject &configuration,
		const std::string &file,
		const std::vector<std::string> &args,
		const std::map<size_t, std::string> &defaultArgs,
		const std::map<std::string, std::string> &variables)
{
	ECFRedParameters redParameters;
	if (environment->MPIrank == 0) {
		std::ifstream ecffile(file);
		redParameters.hadValidECF = ecffile.good();
	}
	int valid = redParameters.hadValidECF;
	MPI_Bcast(&valid, 1, MPI_INT, 0, environment->MPICommunicator);
	redParameters.hadValidECF = valid;

	if (!redParameters.hadValidECF) {
		return redParameters;
	}

	std::vector<std::string> prefix;
	std::vector<std::string> values;
	std::stack<ECFObject*> confStack;
	std::stack<Tokenizer*> tokenStack;

	bool correctlyLoaded = true;
	std::map<size_t, std::vector<std::string> > arguments;

	std::string link;
	ECFParameter *parameter;
	confStack.push(&configuration);
	tokenStack.push(new Tokenizer(file));
	while (tokenStack.size()) {
		switch (tokenStack.top()->next()) {
		case Tokenizer::Token::END:
			delete tokenStack.top();
			tokenStack.pop();
			break;
		case Tokenizer::Token::STRING:
			values.push_back(tokenStack.top()->value());
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
				tokenStack.push(new Tokenizer(value));
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
						ESINFO(GLOBAL_ERROR) << "Invalid argument '" << value << "'";
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
								ESINFO(GLOBAL_ERROR) << "parameter cannot be the [ARG].\n" << tokenStack.top()->lastLines(2);
							}
						}
					}
				}
				break;
			}
			if (value.size() > 4 && StringCompare::caseInsensitiveSuffix(value, ".csv")) {
				Tokenizer csv(value);
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
						ESINFO(GLOBAL_ERROR) << "Error while reading file '" << value << "'";
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
				ESINFO(GLOBAL_ERROR) << "PARSE ERROR: Opening of an unnamed region is not allowed.\n" << tokenStack.top()->lastLines(2);
			}
			if (values.size() > 1) {
				ESINFO(GLOBAL_ERROR) << "PARSE ERROR: Multiple names for a region are not allowed.\n" << tokenStack.top()->lastLines(2);
			}
			prefix.push_back(values[0]);
			parameter = confStack.top()->getParameter(values[0]);
			if (parameter == NULL) {
				ESINFO(GLOBAL_ERROR) << "PARSE ERROR: Unexpected parameter '" << prefix << "'\n" << tokenStack.top()->lastLines(2);
			}
			if (!parameter->isObject()) {
				ESINFO(GLOBAL_ERROR) << "PARSE ERROR: Expected parameter instead of object '" << prefix << "'\n" << tokenStack.top()->lastLines(2);
			}
			redParameters.parameters.push_back(parameter);
			confStack.push(dynamic_cast<ECFObject*>(parameter));
			values.clear();
			break;
		case Tokenizer::Token::OBJECT_CLOSE:
			if (!confStack.size() || !prefix.size()) {
				ESINFO(GLOBAL_ERROR) << "PARSE ERROR: Unexpected region end.\n" << tokenStack.top()->lastLines(2);
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
				ESINFO(GLOBAL_ERROR) << "PARSE ERROR: Incorrect assignment format on line " << tokenStack.top()->line() << ". Use 'PARAMETER' 'VALUE';\n" << tokenStack.top()->lastLines(2);
			}
			std::stringstream ss;
			ss << values[1];
			for (size_t i = 2; i < values.size(); i++) {
				ss << " " << values[i];
			}
			parameter = confStack.top()->getParameter(values[0]);
			if (parameter == NULL) {
				prefix.push_back(values[0]);
				ESINFO(GLOBAL_ERROR) << "PARSE ERROR: Unexpected parameter '" << prefix << "'\n" << tokenStack.top()->lastLines(2);
			}
			if (parameter->isObject()) {
				ESINFO(GLOBAL_ERROR) << "Invalid ECF configuration. Parameter '" << prefix << "::" << parameter->name << "' is an object. Expected '{' instead of '" << ss.str() << "'";
			}
			if (!parameter->setValue(ss.str())) {
				ESINFO(GLOBAL_ERROR) << "PARSE ERROR: Parameter '" << values[0] << "' has wrong value '" << ss.str() << "'";
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
				ESINFO(GLOBAL_ERROR) << "PARSE ERROR: Expected ';' at the end of each expression.\n" << tokenStack.top()->lastLines(1);
			}
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "PARSE ERROR: Unknown token in configuration file";
		}
	}
	if (confStack.size() != 1) {
		ESINFO(GLOBAL_ERROR) << "PARSE ERROR: Unexpected EOF before close all regions.";
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
		ESINFO(GLOBAL_ERROR) << error;
	}
	return redParameters;
}

void ECFReader::set(const Environment &env, const OutputConfiguration &output)
{
	Info::setLevel(env.verbose_level, env.testing_level);
	Measure::setLevel(env.measure_level);
	Logging::path = output.path;
	Logging::debug = env.log_dir;
	Logging::rank = env.MPIrank;
	copyInputData();
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
