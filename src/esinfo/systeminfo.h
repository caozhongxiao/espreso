
#ifndef SRC_ESINFO_SYSTEMINFO_H_
#define SRC_ESINFO_SYSTEMINFO_H_

namespace espreso {
namespace info {
namespace system {

	enum class OPERATIONSYSTEM {
		UNIX
	};

	enum class INSTRUCTIONSET {
		SSE,
		AVX
	};

	enum class BUILD {
		RELEASE,
		MEASUREMENT,
		DEVEL,
		DEBUG
	};

	constexpr OPERATIONSYSTEM os =
		OPERATIONSYSTEM::UNIX;

	constexpr INSTRUCTIONSET instructionSet =
		INSTRUCTIONSET::SSE;

	constexpr BUILD build =
		BUILD::RELEASE;

	void setSignals();
}
}
}


#endif /* SRC_ESINFO_SYSTEMINFO_H_ */
