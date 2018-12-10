
#ifndef SRC_GLOBALS_SYSTEM_H_
#define SRC_GLOBALS_SYSTEM_H_

namespace espreso {

struct system {

	enum class OPERATIONSYSTEM {
		UNIX
	};

	enum class INSTRUCTIONSET {
		SSE,
		AVX
	};

	enum class BUILD {
		DEBUG,
		RELEASE
	};

	static constexpr OPERATIONSYSTEM OS()
	{
		return OPERATIONSYSTEM::UNIX;
	}

	static constexpr INSTRUCTIONSET instructionSet()
	{
		return INSTRUCTIONSET::SSE;
	}

	static constexpr BUILD build()
	{
		return BUILD::DEBUG;
	}

	static void setSignals();
};

}


#endif /* SRC_GLOBALS_SYSTEM_H_ */
