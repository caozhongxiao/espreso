
#ifndef SRC_ESINFO_ESLOG_H_
#define SRC_ESINFO_ESLOG_H_

namespace espreso {
namespace eslog {

class Logger;

void init(const char* ecf, const char* outputPath, int duplication);
void finish();

void start(const char* name, const char* section);
void checkpoint(const char* name);
void end(const char* name);
void ln();
void nextStep(int step);

void startln(const char* name, const char* section);
void checkpointln(const char* name);
void endln(const char* name);

void param(const char* name, const int &value);
void param(const char* name, const long &value);
void param(const char* name, const long unsigned int &value);
void param(const char* name, const double &value);
void param(const char* name, const char* value);

void addsolverparam(const char* name, const char* shortcut, const char* format, int &value);
void addsolverparam(const char* name, const char* shortcut, const char* format, long &value);
void addsolverparam(const char* name, const char* shortcut, const char* format, long unsigned int &value);
void addsolverparam(const char* name, const char* shortcut, const char* format, double &value);
void addsolverparam(const char* name, const char* shortcut, const char* format, const char* &value);
void addsolverparam(const char* name, const char* shortcut, const char* format, bool &value);
void printsolverheader();
void printsolver();

void info(const char* msg);
void solver(const char* msg);
void linearsolver(const char* msg);
void duration(const char* msg);
void warning(const char* msg);
void debug(const char* msg);
void error(const char* msg);
void globalerror(const char* msg);

const char* path();
const char* name();
double time();
double duration();

bool printtime();

extern Logger *logger;

}
}



#endif /* SRC_ESINFO_ESLOG_H_ */
