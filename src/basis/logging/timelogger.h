
#ifndef SRC_BASIS_LOGGING_TIMELOGGER_H_
#define SRC_BASIS_LOGGING_TIMELOGGER_H_

#include "verbosity.h"

#include <vector>

namespace espreso {

class TimeLogger: public Verbosity<TimeLogger, 't'> {

	struct Event {
		const char* name;

		union Data {
			double time;
			int    ivalue;
			long   lvalue;
			size_t svalue;
			double dvalue;
		} data;

		enum {
			START, CHECKPOINT, END,
			INT, LONG, SIZE, DOUBLE
		} type;
	};

public:
	static double time();

	void start(const char* region)
	{
		_events.push_back(Event{ region, Event::Data{ .time = time() }, Event::START });
	}

	void checkpoint(const char* region)
	{
		_events.push_back(Event{ region, Event::Data{ .time = time() }, Event::CHECKPOINT });
	}

	void end(const char* region)
	{
		_events.push_back(Event{ region, Event::Data{ .time = time() }, Event::END });
	}

	void param(const char* name, const int &value)
	{
		_events.push_back(Event{ name, Event::Data{ .ivalue = value }, Event::INT });
	}

	void param(const char* name, const long &value)
	{
		_events.push_back(Event{ name, Event::Data{ .lvalue = value }, Event::LONG });
	}

	void param(const char* name, const long unsigned int &value)
	{
		_events.push_back(Event{ name, Event::Data{ .svalue = value }, Event::SIZE });
	}

	void param(const char* name, const double &value)
	{
		_events.push_back(Event{ name, Event::Data{ .dvalue = value }, Event::DOUBLE });
	}

	void param(const char* name, const char* value)
	{
		// do nothing
	}

	void ln()
	{
		// do nothing
	}

	TimeLogger()
	{
		_init = time();
		_events.reserve(1000000);
	}

	void evaluate();

protected:
	double _init;
	std::vector<Event> _events;
};

}



#endif /* SRC_BASIS_LOGGING_TIMELOGGER_H_ */
