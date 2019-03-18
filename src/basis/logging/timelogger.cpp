
#include "timelogger.h"
#include "progresslogger.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/systeminfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/timeinfo.h"
#include "basis/utilities/parser.h"
#include "basis/utilities/communication.h"
#include "input/input.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/elementsregionstore.h"
#include "mesh/store/boundaryregionstore.h"

#include "omp.h"
#include <cstdio>
#include <cstring>

using namespace espreso;

double TimeLogger::init = TimeLogger::time();

double TimeLogger::time()
{
	return omp_get_wtime();
}

double TimeLogger::duration()
{
	return omp_get_wtime() - init;
}

void TimeLogger::mergeEvents(void *in, void *out, int *len, MPI_Datatype *datatype)
{
	int size = *len / sizeof(TimeLogger::EventStatistics);
	for (int i = 0; i < size; i++) {
		Event::Data &inmin = (static_cast<TimeLogger::EventStatistics*>(in) + i)->min;
		Event::Data &inmax = (static_cast<TimeLogger::EventStatistics*>(in) + i)->max;
		Event::Data &inavg = (static_cast<TimeLogger::EventStatistics*>(in) + i)->avg;
		Event::Data &indmin = (static_cast<TimeLogger::EventStatistics*>(in) + i)->dmin;
		Event::Data &indmax = (static_cast<TimeLogger::EventStatistics*>(in) + i)->dmax;
		Event::Data &indavg = (static_cast<TimeLogger::EventStatistics*>(in) + i)->davg;
		Event::Data &outmin = (static_cast<TimeLogger::EventStatistics*>(out) + i)->min;
		Event::Data &outmax = (static_cast<TimeLogger::EventStatistics*>(out) + i)->max;
		Event::Data &outavg = (static_cast<TimeLogger::EventStatistics*>(out) + i)->avg;
		Event::Data &outdmin = (static_cast<TimeLogger::EventStatistics*>(out) + i)->dmin;
		Event::Data &outdmax = (static_cast<TimeLogger::EventStatistics*>(out) + i)->dmax;
		Event::Data &outdavg = (static_cast<TimeLogger::EventStatistics*>(out) + i)->davg;

		switch ((static_cast<TimeLogger::EventStatistics*>(out) + i)->type) {
		case Event::START:
		case Event::CHECKPOINT:
		case Event::END:
			outmin.time = std::min(outmin.time, inmin.time);
			outmax.time = std::max(outmax.time, inmax.time);
			outavg.time += inavg.time;
			outdmin.time = std::min(outdmin.time, indmin.time);
			outdmax.time = std::max(outdmax.time, indmax.time);
			outdavg.time += indavg.time;
			break;
		case Event::INT:
			outmin.ivalue = std::min(outmin.ivalue, inmin.ivalue);
			outmax.ivalue = std::max(outmax.ivalue, inmax.ivalue);
			outavg.ivalue += inavg.ivalue;
			break;
		case Event::LONG:
			outmin.lvalue = std::min(outmin.lvalue, inmin.lvalue);
			outmax.lvalue = std::max(outmax.lvalue, inmax.lvalue);
			outavg.lvalue += inavg.lvalue;
			break;
		case Event::SIZE:
			outmin.svalue = std::min(outmin.svalue, inmin.svalue);
			outmax.svalue = std::max(outmax.svalue, inmax.svalue);
			outavg.svalue += inavg.svalue;
			break;
		case Event::DOUBLE:
			outmin.dvalue = std::min(outmin.dvalue, inmin.dvalue);
			outmax.dvalue = std::max(outmax.dvalue, inmax.dvalue);
			outavg.dvalue += inavg.dvalue;
			break;
		default:
			break;
		}
	}
}

static void printdata(
		ProgressLogger &logger,
		const char* format, const char* name, const char* suffix,
		double avg, double min, double max, double sectiontime)
{
	std::string fullname = std::string(name) + std::string(suffix);
	avg /= info::mpi::size;
	std::string savg = std::to_string(avg);
	std::string smin = std::to_string(min);
	std::string smax = std::to_string(max);
	savg[8] = smin[8] = smax[8] = '\0';
	logger.info(format, fullname.c_str(), savg.c_str(), smin.c_str(), smax.c_str(), 100 * avg / sectiontime, max / min);
}

void TimeLogger::evaluate(ProgressLogger &logger)
{
	std::vector<EventStatistics> events(_events.begin(), _events.end());
	std::vector<double> prev(10), begin(10);

	size_t namewidth = 43, width = 78;

	for (size_t i = 0; i < events.size(); i++) {
		switch (events[i].type) {
		case Event::START:
			prev.push_back(events[i].data.time);
			begin.push_back(events[i].data.time);
			events[i].data.time -= init;
			events[i].duration.time = 0;
			break;
		case Event::CHECKPOINT:
			events[i].data.time -= prev.back();
			events[i].duration.time -= begin.back();
			prev.back() += events[i].data.time;
			namewidth = std::max(namewidth, strlen(events[i].name));
			break;
		case Event::END:
			events[i].data.time -= prev.back();
			events[i].duration.time -= begin.back();
			prev.pop_back();
			begin.pop_back();
			namewidth = std::max(namewidth, strlen(events[i].name));
			break;
		default:
			namewidth = std::max(namewidth, strlen(events[i].name) + 2);
			break;
		}
		events[i].min = events[i].data;
		events[i].max = events[i].data;
		events[i].avg = events[i].data;
		events[i].dmin = events[i].duration;
		events[i].dmax = events[i].duration;
		events[i].davg = events[i].duration;
	}

	std::vector<EventStatistics> statistics(events);

	{ // synchroniza across processes
		size_t eventsize = events.size(), minsize, maxsize;
		MPI_Allreduce(&eventsize, &minsize, sizeof(size_t), MPI_BYTE, MPITools::operations->SIZET.min, info::mpi::comm);
		MPI_Allreduce(&eventsize, &maxsize, sizeof(size_t), MPI_BYTE, MPITools::operations->SIZET.max, info::mpi::comm);
		if (minsize == eventsize && maxsize == eventsize) {
			MPI_Op reduce;
			MPI_Op_create(&TimeLogger::mergeEvents, 1, &reduce);
			MPI_Reduce(events.data(), statistics.data(), sizeof(EventStatistics) * events.size(), MPI_BYTE, reduce, 0, info::mpi::comm);
			MPI_Op_free(&reduce);
		} else {
			logger.warning("Various number of time events (only root data are printed)\n");
		}
	}

	if (info::mpi::rank) {
		MPI_Barrier(info::mpi::comm);
		return;
	}

	// print coarse run and mesh statistics
	auto totalesize = [] (std::vector<esint> &ecounters) {
		esint size = 0;
		for (size_t etype = 0; etype < ecounters.size(); etype++) {
			size += ecounters[etype];
		}
		return size;
	};

	size_t eregs = info::mesh->elementsRegions.size() - 1;
	size_t bregs = 0;
	size_t nregs = 0;
	if (StringCompare::caseInsensitiveEq(info::mesh->elementsRegions.back()->name, "NAMELESS_ELEMENT_SET")) {
		--eregs;
	}
	for (size_t i = 1; i < info::mesh->boundaryRegions.size(); i++) {
		if (info::mesh->boundaryRegions[i]->dimension) {
			++bregs;
		} else {
			++nregs;
		}
	}

	logger.info(" ============================================================================================= \n");
	logger.info(" == commit   %*s == \n", width, info::system::commit());
	logger.info(" == ecf      %*s == \n", width, info::ecf->ecffile.c_str());
	logger.info(" == MPI      %*d == \n", width, info::mpi::size);
	logger.info(" == OMP      %*d == \n", width, info::env::OMP_NUM_THREADS);
	logger.info(" == mesh     %*s == \n", width, Input::inputFile(*info::ecf));
	logger.info(" == elements %*d == \n", width, totalesize(info::mesh->elements->ecounters));
	logger.info(" == nodes    %*d == \n", width, info::mesh->nodes->uniqueTotalSize);
	logger.info(" == elements regions %*ld == \n", width - 8, eregs);
	logger.info(" == boundary regions %*ld == \n", width - 8, bregs);
	logger.info(" == nodes regions    %*ld == \n", width - 8, nregs);
	logger.info(" ============================================================================================= \n");

	double duration = TimeLogger::duration();
	const char* headformat = " %-44s %s  <%s - %s> [%5.2f] [%5.2f]\n";
	const char* dataformat = "  %-43s %s  <%s - %s> [%5.2f] [%5.2f]\n";

	auto print = [&] (size_t start, size_t end, int printeddepth, const std::vector<const char*> &duplicities) {
		logger.info(" ============================================ avg. [s]  < min [s] -  max [s]> [  %%  ] [ imb ]   \n");
		int depth = printeddepth - 1;
		std::vector<const char*> printed;
		for (size_t i = start; i <= end; i++) {
			switch (statistics[i].type) {
			case Event::START:
				++depth;
				if (depth == printeddepth) {
					printdata(logger, headformat, events[i].name, " START AT",
							statistics[i].avg.time, statistics[i].min.time, statistics[i].max.time, duration);
				}
				break;
			case Event::CHECKPOINT:
				if (depth == printeddepth) {
					auto it = std::find(duplicities.begin(), duplicities.end(), events[i].name);
					if (it == duplicities.end()) {
						printdata(logger, dataformat, events[i].name, "",
								statistics[i].avg.time, statistics[i].min.time, statistics[i].max.time, statistics[end].duration.time);
					} else {
						if (std::find(printed.begin(), printed.end(), events[i].name) == printed.end()) {
							printed.push_back(events[i].name);
							int count = 0;
							double avg = 0, min = duration, max = 0;
							for (size_t j = start; j <= end; j++) {
								if (events[i].name == events[j].name) {
									if (count) {
										min = std::min(min, statistics[j].min.time);
										max = std::max(max, statistics[j].max.time);
										avg += statistics[j].avg.time;
									}
									++count;
								}
							}
							logger.info("  %s [%dx]\n", events[i].name, count);
							printdata(logger, dataformat, "  [run=FIRST]", "", statistics[i].avg.time, statistics[i].min.time, statistics[i].max.time, statistics[end].duration.time);
							printdata(logger, dataformat, "  [run=REST]", "", avg / (count - 1), min, max, statistics[end].duration.time);
						}
					}
				}
				break;
			case Event::END:
				if (depth == printeddepth) {
					printdata(logger, dataformat, events[i].name, "",
							statistics[i].avg.time, statistics[i].min.time, statistics[i].max.time, statistics[end].duration.time);
					printdata(logger, headformat, events[start].name, " TOTAL DURATION",
							statistics[i].davg.time, statistics[i].dmin.time, statistics[i].dmax.time, duration);
				}
				--depth;
				break;
			case Event::INT:
				if (depth > 1 && depth == printeddepth) {
					logger.info("    [param=%s] %*d  <%8d - %8d>         [%5.2f]\n",
							events[i].name, namewidth - strlen(events[i].name) - 2,
							statistics[i].avg.ivalue / info::mpi::size,
							statistics[i].min.ivalue, statistics[i].max.ivalue,
							(double)statistics[i].max.ivalue / statistics[i].min.ivalue);
				}
				break;
			case Event::LONG:
				if (depth > 1 && depth == printeddepth) {
					logger.info("    [param=%s] %*ld  <%8ld - %8ld>         [%5.2f]\n",
							events[i].name, namewidth - strlen(events[i].name) - 2,
							statistics[i].avg.lvalue / info::mpi::size,
							statistics[i].min.lvalue, statistics[i].max.lvalue,
							(double)statistics[i].max.lvalue / statistics[i].min.lvalue);
				}
				break;
			case Event::SIZE:
				if (depth > 1 && depth == printeddepth) {
					logger.info("    [param=%s] %*ld  <%8ld - %8ld>         [%5.2f]\n",
							events[i].name, namewidth - strlen(events[i].name) - 2,
							statistics[i].avg.svalue / info::mpi::size,
							statistics[i].min.svalue, statistics[i].max.svalue,
							(double)statistics[i].max.svalue / statistics[i].min.svalue);
				}
				break;
			case Event::DOUBLE:
				if (depth > 1 && depth == printeddepth) {
					logger.info("    [param=%s] %*f  <%*f - %*f>         [%5.2f]\n",
							events[i].name, namewidth - strlen(events[i].name) - 2,
							statistics[i].avg.dvalue / info::mpi::size, 8,
							statistics[i].min.dvalue, 8, statistics[i].max.dvalue,
							statistics[i].max.dvalue / statistics[i].min.dvalue);
				}
				break;
			default:
				break;
			}
		}
		logger.info(" ============================================================================================= \n");
	};

	// WARNING: MPI sometimes copy name pointer from other process, hence use _events names
	std::vector<size_t> begins;
	std::vector<std::vector<const char*> > uniques, duplicities;
	std::vector<std::string> solvers(time::step + 1);
	int loadstep = 0;
	size_t lastend = 0;

	auto addloadstep = [&] (size_t i) {
		if (begins.size() == 1 && strcmp(events[i].name, "ESPRESO: SOLVED") == 0) {
			if (loadstep > 9) {
				solvers[loadstep] = "ESPRESO: SOLVED [LOADSTEP " + std::to_string(loadstep) + "]";
			} else {
				solvers[loadstep] = "ESPRESO: SOLVED [LOADSTEP  " + std::to_string(loadstep) + "]";
			}
			events[i].name = solvers[loadstep].c_str();
		}
	};

	for (size_t i = 0; i < statistics.size(); i++) {
		switch (statistics[i].type) {
		case Event::START:
			begins.push_back(i);
			uniques.push_back({});
			duplicities.push_back({});
			break;
		case Event::CHECKPOINT:
			addloadstep(i);
			if (std::find(uniques.back().begin(), uniques.back().end(), events[i].name) == uniques.back().end()) {
				uniques.back().push_back(events[i].name);
			} else {
				if (std::find(duplicities.back().begin(), duplicities.back().end(), events[i].name) == duplicities.back().end()) {
					duplicities.back().push_back(events[i].name);
				}
			}
			break;
		case Event::END:
			if (begins.size() > 1) {
				print(begins.back(), i, begins.size() + 1, duplicities.back()); // do not print loadsteps
			}
			addloadstep(i);
			begins.pop_back();
			uniques.pop_back();
			duplicities.pop_back();
			lastend = i;
			break;
		case Event::LOADSTEP:
			++loadstep;
			if (verbosity > 1) {
				logger.info(" == LOADSTEP: %2d ============================================================================= \n", loadstep);
			}
			break;
		default:
			break;
		}
	}

	logger.info(" == OVERALL TIME ============================================================================= \n");
	print(0, lastend, 0, {});
	MPI_Barrier(info::mpi::comm);
}



