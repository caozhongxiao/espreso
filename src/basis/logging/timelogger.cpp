
#include "timelogger.h"
#include "progresslogger.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/systeminfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
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

double TimeLogger::time()
{
	return omp_get_wtime();
}

void TimeLogger::mergeEvents(void *in, void *out, int *len, MPI_Datatype *datatype)
{
	int size = *len / sizeof(TimeLogger::EventStatistics);
	for (int i = 0; i < size; i++) {
		Event::Data &inmin = (static_cast<TimeLogger::EventStatistics*>(in) + i)->min;
		Event::Data &inmax = (static_cast<TimeLogger::EventStatistics*>(in) + i)->max;
		Event::Data &inavg = (static_cast<TimeLogger::EventStatistics*>(in) + i)->avg;
		Event::Data &outmin = (static_cast<TimeLogger::EventStatistics*>(out) + i)->min;
		Event::Data &outmax = (static_cast<TimeLogger::EventStatistics*>(out) + i)->max;
		Event::Data &outavg = (static_cast<TimeLogger::EventStatistics*>(out) + i)->avg;

		switch ((static_cast<TimeLogger::EventStatistics*>(out) + i)->type) {
		case Event::START:
		case Event::CHECKPOINT:
		case Event::END:
			outmin.time = std::min(outmin.time, inmin.time);
			outmax.time = std::max(outmax.time, inmax.time);
			outavg.time += inavg.time;
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

void TimeLogger::evaluate(ProgressLogger &logger)
{
	double duration = time() - _init;
	std::vector<double> prev(10);

	size_t namewidth = 44, width = 78;

	for (size_t i = 0; i < _events.size(); i++) {
		switch (_events[i].type) {
		case Event::START:
			prev.push_back(_events[i].data.time);
			_events[i].data.time -= _init;
			break;
		case Event::CHECKPOINT:
			_events[i].data.time -= prev.back();
			prev.back() += _events[i].data.time;
			namewidth = std::max(namewidth, strlen(_events[i].name));
			break;
		case Event::END:
			_events[i].data.time -= prev.back();
			prev.pop_back();
			namewidth = std::max(namewidth, strlen(_events[i].name));
			break;
		default:
			namewidth = std::max(namewidth, strlen(_events[i].name) + 2);
			break;
		}
	}
	std::vector<EventStatistics> events(_events.begin(), _events.end());
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
			if (info::mpi::rank == 0) {
				logger.warning("Various number of time events (only root data are printed)\n");
			}
		}
	}

	if (info::mpi::rank) {
		MPI_Barrier(info::mpi::comm);
		return;
	}

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

	auto print = [&] (size_t start, size_t end, int printeddepth, int loadstep) {
		logger.info(" ============================================ avg. [s]  < min [s] -  max [s]> [  %%  ] [ imb ]   \n");
		int depth = printeddepth - 1;
		for (size_t i = start; i < end; i++) {
			switch (statistics[i].type) {
			case Event::START:
				++depth;
				--namewidth;
				if (depth == printeddepth) {
					logger.info("  %-*s %f  <%f - %f> [-----] [%5.2f]\n",
							namewidth, _events[i].name,
							statistics[i].avg.time / info::mpi::size,
							statistics[i].min.time, statistics[i].max.time,
							statistics[i].max.time / statistics[i].min.time);
				}
				break;
			case Event::CHECKPOINT:
				if (depth == printeddepth) {
					if (loadstep > 0) {
						logger.info("  %s [LOADSTEP:%2d] %*f  <%f - %f> [%5.2f] [%5.2f]\n",
									_events[i].name, loadstep, namewidth - 21,
									statistics[i].avg.time / info::mpi::size,
									statistics[i].min.time, statistics[i].max.time,
									100 * (statistics[i].avg.time / info::mpi::size) / duration,
									statistics[i].max.time / statistics[i].min.time);
					} else {
						logger.info("  %-*s %f  <%f - %f> [%5.2f] [%5.2f]\n",
									namewidth, _events[i].name,
									statistics[i].avg.time / info::mpi::size,
									statistics[i].min.time, statistics[i].max.time,
									100 * (statistics[i].avg.time / info::mpi::size) / duration,
									statistics[i].max.time / statistics[i].min.time);
					}
				}
				break;
			case Event::END:
				if (depth == printeddepth) {
					if (loadstep > 0) {
						logger.info("  %s [LOADSTEP:%2d] %*f  <%f - %f> [%5.2f] [%5.2f]\n",
									_events[i].name, loadstep, namewidth - 21,
									statistics[i].avg.time / info::mpi::size,
									statistics[i].min.time, statistics[i].max.time,
									100 * (statistics[i].avg.time / info::mpi::size) / duration,
									statistics[i].max.time / statistics[i].min.time);
					} else {
						logger.info("  %-*s %f  <%f - %f> [%5.2f] [%5.2f]\n",
									namewidth, _events[i].name,
									statistics[i].avg.time / info::mpi::size,
									statistics[i].min.time, statistics[i].max.time,
									100 * (statistics[i].avg.time / info::mpi::size) / duration,
									statistics[i].max.time / statistics[i].min.time);
					}
				}
				--depth;
				++namewidth;
				break;
			case Event::LOADSTEP:
				++loadstep;
				break;
			case Event::INT:
				if (depth > 1 && depth == printeddepth) {
					logger.info("    [param=%s] %*d  <%8d - %8d>         [%5.2f]\n",
							_events[i].name, namewidth - strlen(_events[i].name) - 2,
							statistics[i].avg.ivalue / info::mpi::size,
							statistics[i].min.ivalue, statistics[i].max.ivalue,
							(double)statistics[i].max.ivalue / statistics[i].min.ivalue);
				}
				break;
			case Event::LONG:
				if (depth > 1 && depth == printeddepth) {
					logger.info("    [param=%s] %*ld  <%8ld - %8ld>         [%5.2f]\n",
							_events[i].name, namewidth - strlen(_events[i].name) - 2,
							statistics[i].avg.lvalue / info::mpi::size,
							statistics[i].min.lvalue, statistics[i].max.lvalue,
							(double)statistics[i].max.lvalue / statistics[i].min.lvalue);
				}
				break;
			case Event::SIZE:
				if (depth > 1 && depth == printeddepth) {
					logger.info("    [param=%s] %*ld  <%8ld - %8ld>         [%5.2f]\n",
							_events[i].name, namewidth - strlen(_events[i].name) - 2,
							statistics[i].avg.svalue / info::mpi::size,
							statistics[i].min.svalue, statistics[i].max.svalue,
							(double)statistics[i].max.svalue / statistics[i].min.svalue);
				}
				break;
			case Event::DOUBLE:
				if (depth > 1 && depth == printeddepth) {
					logger.info("    [param=%s] %*f  <%*f - %*f>         [%5.2f]\n",
							_events[i].name, namewidth - strlen(_events[i].name) - 2,
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

	// WARNING: MPI sometimes copy name pointer from other process, hence use _events names
	std::vector<size_t> begins;
	int loadstep = 0;
	for (size_t i = 0; i < statistics.size(); i++) {
		switch (statistics[i].type) {
		case Event::START:
			begins.push_back(i);
			break;
		case Event::END:
			if (begins.size() > 1) {
				print(begins.back(), i + 1, begins.size() + 1, -100); // do not print loadsteps
			}
			begins.pop_back();
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
	print(0, statistics.size(), 0, 0);
	MPI_Barrier(info::mpi::comm);
}



