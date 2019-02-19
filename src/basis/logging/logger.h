
#ifndef SRC_BASIS_LOGGING_LOGGER_H_
#define SRC_BASIS_LOGGING_LOGGER_H_

#include <cstdio>
#include <type_traits>

namespace espreso {

#define __ES__PACK(...) __VA_ARGS__
#define __ES__FORALL(fnc, header, call)        \
	template<class Logger, class... Other>     \
	typename                                   \
	std::enable_if<sizeof...(Other)>::type     \
	fnc(header)                                \
	{                                          \
		fnc<Logger>(call);                     \
		fnc<Other...>(call);                   \
	}                                          \
	void fnc(header)                           \
	{                                          \
		fnc<Loggers...>(call);                 \
	}

template <class... Loggers>
class Logger: public Loggers... {

public:
	template<class Logger>
	void start(const char* n)
	{
		++Logger::level;
		if (Logger::isAllowed()) {
			Logger::start(n);
		}
	}

	template<class Logger>
	void checkpoint(const char* n)
	{
		if (Logger::isAllowed()) {
			Logger::checkpoint(n);
		}
	}

	template<class Logger>
	void end(const char* n)
	{
		if (Logger::isAllowed()) {
			Logger::end(n);
		}
		Logger::finishing = true;
	}

	template<class Logger>
	void param(const char* n, const int &value)
	{
		if (Logger::isAllowed()) {
			Logger::param(n, value);
		}
	}

	template<class Logger>
	void param(const char* n, const long &value)
	{
		if (Logger::isAllowed()) {
			Logger::param(n, value);
		}
	}

	template<class Logger>
	void param(const char* n, const long unsigned int &value)
	{
		if (Logger::isAllowed()) {
			Logger::param(n, value);
		}
	}

	template<class Logger>
	void param(const char* n, const double &value)
	{
		if (Logger::isAllowed()) {
			Logger::param(n, value);
		}
	}

	template<class Logger>
	void param(const char* n, const char* value)
	{
		if (Logger::isAllowed()) {
			Logger::param(n, value);
		}
	}

	template<class Logger>
	void ln()
	{
		if (Logger::isAllowed()) {
			Logger::ln();
		}
		if (Logger::finishing) {
			Logger::finishing = false;
			--Logger::level;
		}
	}

	__ES__FORALL(start, const char* n, n);
	__ES__FORALL(checkpoint, const char* n, n);
	__ES__FORALL(end, const char* n, n);
	__ES__FORALL(ln, , );
	__ES__FORALL(param, __ES__PACK(const char* n, const int &value), __ES__PACK(n, value));
	__ES__FORALL(param, __ES__PACK(const char* n, const long &value), __ES__PACK(n, value));
	__ES__FORALL(param, __ES__PACK(const char* n, const long unsigned int &value), __ES__PACK(n, value));
	__ES__FORALL(param, __ES__PACK(const char* n, const double &value), __ES__PACK(n, value));
	__ES__FORALL(param, __ES__PACK(const char* n, const char* value), __ES__PACK(n, value));
};
}



#endif /* SRC_BASIS_LOGGING_LOGGER_H_ */
