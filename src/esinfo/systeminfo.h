
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

	constexpr OPERATIONSYSTEM os = OPERATIONSYSTEM::UNIX;

	constexpr INSTRUCTIONSET instructionSet = INSTRUCTIONSET::SSE;

#ifdef __ESMODE__
	constexpr BUILD build = BUILD::__ESMODE__;
#else
	constexpr BUILD build = BUILD::RELEASE;
#endif


	extern const char* commit;

	void setSignals();
}
}
}


#endif /* SRC_ESINFO_SYSTEMINFO_H_ */
