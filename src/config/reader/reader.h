
#ifndef SRC_CONFIG_READER_READER_H_
#define SRC_CONFIG_READER_READER_H_

#include <string>
#include <vector>
#include <map>

namespace espreso {

struct ECFObject;
struct Environment;
struct OutputConfiguration;

class Reader {

public:
	static void read(
			ECFObject &configuration,
			const std::string &file,
			const std::map<size_t, std::string> &defaultArgs = {},
			const std::map<std::string, std::string> &variables = {}) { _read(configuration, file, {}, defaultArgs, variables); }

	static void read(
			ECFObject &configuration,
			int* argc,
			char ***argv,
			const std::map<size_t, std::string> &defaultArgs = {},
			const std::map<std::string, std::string> &variables = {}) { _read(configuration, argc, argv, defaultArgs, variables); }

	static void set(const Environment &env, const OutputConfiguration &output);

	static void print(const ECFObject &configuration);
	static void store(const ECFObject &configuration, const std::vector<std::string> &subConfigurations);
	static void copyInputData();

private:
	static std::string configurationFile;

	static void _read(
			ECFObject &configuration,
			const std::string &file,
			const std::vector<std::string> &args,
			const std::map<size_t, std::string> &defaultArgs,
			const std::map<std::string, std::string> &variables);

	static void _read(
			ECFObject &configuration,
			int* argc,
			char ***argv,
			const std::map<size_t, std::string> &defaultArgs,
			const std::map<std::string, std::string> &variables);
};

}

#endif /* SRC_CONFIG_READER_READER_H_ */
