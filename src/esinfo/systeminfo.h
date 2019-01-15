
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
		DEBUG,
		MEASUREMENT
	};

	constexpr OPERATIONSYSTEM os =
		OPERATIONSYSTEM::UNIX;

	constexpr INSTRUCTIONSET instructionSet =
		INSTRUCTIONSET::SSE;

	constexpr BUILD build =
		BUILD::DEBUG;

	void setSignals();
}
}
}


#endif /* SRC_ESINFO_SYSTEMINFO_H_ */
