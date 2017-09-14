
#ifndef SRC_CONFIG_READER_READER_H_
#define SRC_CONFIG_READER_READER_H_

#include <string>
#include <vector>
#include <map>
#include <ostream>

namespace espreso {

struct ECFParameter;
struct ECFObject;
struct Environment;
struct OutputConfiguration;

struct ECFRedParameters {
	std::vector<ECFParameter*> parameters;
	std::map<const ECFParameter*, std::string> defaulted;
};

class ECFReader {

public:
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

	static void set(const Environment &env, const OutputConfiguration &output);

	static void store(const ECFObject &configuration, std::ostream &os, bool onlyAllowed = true, bool printPatterns = false, const ECFRedParameters &parameters = {});

	static std::string configurationFile;

private:
	static void copyInputData();
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
