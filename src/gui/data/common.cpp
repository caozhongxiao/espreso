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

std::vector<std::string> Common::fnVariablesParts()
{
    std::vector<std::string> v = fnVariables();
    v.push_back("t");
    v.push_back("ti");
    v.push_back("te");
    v.push_back("tim");
    v.push_back("tem");
    v.push_back("temp");
    v.push_back("tempe");
    v.push_back("temper");
    v.push_back("tempera");
    v.push_back("tempera");
    v.push_back("temperat");
    v.push_back("temperatu");
    v.push_back("temperatur");

    return v;
}
