
#ifndef SRC_CONFIG_READER_READER_H_
#define SRC_CONFIG_READER_READER_H_

#include <string>
#include <vector>
#include <map>
#include <ostream>

namespace espreso {

struct ECFParameter;
struct ECFObject;
struct EnvironmentConfiguration;
struct OutputConfiguration;
struct VerboseArg;

struct ECFRedParameters {
	bool hadValidECF;
	std::vector<ECFParameter*> parameters;
	std::map<const ECFParameter*, std::string> defaulted;
};

class ECFReader {

public:
	ECFReader(int argc, char **argv);
	ECFReader(const std::string &file);

	std::string ecf;
	std::string outputPath;
	int meshDuplication;

	static ECFRedParameters read(
			ECFObject &configuration,
			const std::string &file,
			const std::map<size_t, std::string> &defaultArgs = {},
			const std::map<std::string, std::string> &variables = {})
	{
		return _read(configuration, file, {}, defaultArgs, variables);
	}

	static ECFRedParameters read(
			ECFObject &configuration,
			int* argc,
			char ***argv,
			const std::map<size_t, std::string> &defaultArgs = {},
			const std::map<std::string, std::string> &variables = {})
	{
		return _read(configuration, argc, argv, defaultArgs, variables);
	}

	static void store(const ECFObject &configuration, std::ostream &os, bool onlyAllowed = true, bool printPatterns = false, const ECFRedParameters &parameters = ECFRedParameters());

	static std::string configurationFile;
	static std::vector<VerboseArg*> verbosity;

private:
	static ECFRedParameters _read(
			ECFObject &configuration,
			const std::string &file,
			const std::vector<std::string> &args,
			const std::map<size_t, std::string> &defaultArgs,
			const std::map<std::string, std::string> &variables);

	static ECFRedParameters _read(
			ECFObject &configuration,
			int* argc,
			char ***argv,
			const std::map<size_t, std::string> &defaultArgs,
			const std::map<std::string, std::string> &variables);
};

}

#endif /* SRC_CONFIG_READER_READER_H_ */
