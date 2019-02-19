
#ifndef SRC_BASIS_LOGGING_VERBOSITY_H_
#define SRC_BASIS_LOGGING_VERBOSITY_H_

namespace espreso {

struct VerboseArg {
	char argflag;
	int level;
	int verbosity;
	int finishing;

	VerboseArg(char argflag);

	bool isAllowed() const {
		return level <= verbosity;
	}
};

template <typename T, char C>
struct Verbosity: public VerboseArg {
	Verbosity(): VerboseArg(C) {}
};

}



#endif /* SRC_BASIS_LOGGING_VERBOSITY_H_ */
