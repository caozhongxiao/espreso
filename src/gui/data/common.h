#ifndef COMMON_H
#define COMMON_H

#include <string>
#include <vector>

#define REGEXPR_DOUBLE "^[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?$"
#define REGEXPR_VARNAME "^[a-zA-Z_][a-zA-Z_0-9]*$"
#define REGEXPR_NAME "^[A-Za-z]+(\\s+[A-Za-z]+)*(\\s+[0-9]+)?$"

class Common
{
private:
    Common() {}

public:
    static std::vector<std::string> fnVariables();
};

#endif // COMMON_H
