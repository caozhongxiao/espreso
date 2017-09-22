#ifndef COMMON_H
#define COMMON_H

#include <string>
#include <vector>

#define GUI_REGEXPR_DOUBLE "^[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?$"
#define GUI_REGEXPR_VARNAME "^[a-zA-Z_][a-zA-Z_0-9]*$"
#define GUI_REGEXPR_NAME "^[A-Za-z]+(\\s+[A-Za-z]+)*(\\s+[0-9]+)?$"

class Common
{
private:
    Common() {}

public:
    static std::vector<std::string> fnVariables();
    static std::vector<std::string> fnVariablesParts();
};

#endif // COMMON_H
