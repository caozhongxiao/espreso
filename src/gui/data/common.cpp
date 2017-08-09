#include "common.h"

std::vector<std::string> Common::fnVariables()
{
    std::vector<std::string> v;
    v.push_back("x");
    v.push_back("y");
    v.push_back("z");
    v.push_back("time");
    v.push_back("temperature");

    return v;
}
