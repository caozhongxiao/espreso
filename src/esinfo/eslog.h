
#ifndef SRC_ESINFO_ESLOG_H_
#define SRC_ESINFO_ESLOG_H_

namespace espreso {
namespace eslog {

class Logger;

// the constructor of the logger
// it allows to compute time to initialization
void create();
// it allows calling of the rest methods
// we have to get output directory from configuration file!!!
void init(int *argc, char ***argv);
// print the overall statistics and destroy the logger
void finish();

void start(const char* name);
void checkpoint(const char* name);
void end(const char* name);
void ln();

void startln(const char* name);
void checkpointln(const char* name);
void endln(const char* name);

void param(const char* name, const int &value);
void param(const char* name, const long &value);
void param(const char* name, const long unsigned int &value);
void param(const char* name, const double &value);
void param(const char* name, const char* value);

const char* path();
const char* name();

bool printtime();

extern Logger *logger;

}
}



#endif /* SRC_ESINFO_ESLOG_H_ */
