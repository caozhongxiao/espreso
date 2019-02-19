
#include "verbosity.h"

#include "config/reader/reader.h"

using namespace espreso;

VerboseArg::VerboseArg(char argflag)
: argflag(argflag), level(0), verbosity(1), finishing(0)
{
	ECFReader::verbosity.push_back(this);
}


