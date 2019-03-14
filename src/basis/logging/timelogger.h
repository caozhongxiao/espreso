
#ifndef SRC_BASIS_LOGGING_TIMELOGGER_H_
#define SRC_BASIS_LOGGING_TIMELOGGER_H_

#include "verbosity.h"

#include "mpi.h"
#include <vector>


namespace espreso {

class ProgressLogger;

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
			START, CHECKPOINT, END, LOADSTEP,
			INT, LONG, SIZE, DOUBLE
		} type;
	};

	struct EventStatistics: public Event {
		Data min, max, avg;
		Data duration, dmin, dmax, davg;

		EventStatistics(const Event &event)
		: Event(event),
		  min(event.data), max(event.data), avg(event.data),
		  duration(event.data), dmin(event.data), dmax(event.data), davg(event.data)
		{ }
	};

public:
	static double time();
	static double duration();

	void start(const char* region, const char* section)
	{
		_events.push_back(Event{ section, Event::Data{ .time = time() }, Event::START });
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

	void nextLoadStep(int step)
	{
		_events.push_back(Event{ "NEXT STEP", Event::Data{ .ivalue = step }, Event::LOADSTEP });
	}

	TimeLogger()
	{
		_events.reserve(1000000);
	}

	void evaluate(ProgressLogger &logger);

protected:
	static void mergeEvents(void *in, void *out, int *len, MPI_Datatype *datatype);

	static double init;
	std::vector<Event> _events;
};

}



#endif /* SRC_BASIS_LOGGING_TIMELOGGER_H_ */
