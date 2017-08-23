
#ifndef SRC_CONFIG_VALUEHOLDER_H_
#define SRC_CONFIG_VALUEHOLDER_H_

#include "configuration.h"
#include "../basis/utilities/parser.h"
#include <sstream>

namespace espreso {

template <typename Ttype>
struct ECFValueHolder: public ECFValue {
	Ttype &value;

	ECFValueHolder(Ttype &value): value(value) {}

	std::string getValue() const
	{
		std::stringstream ss;
		ss << value;
		return ss.str();
	}

	bool setValue(const std::string &value)
	{
		std::stringstream ss(value);
		ss >> this->value;
		return ss.eof() && !ss.fail();
	}
};

template <>
inline std::string ECFValueHolder<std::string>::getValue() const
{
	return value;
}

template <>
inline std::string ECFValueHolder<ECFExpression>::getValue() const
{
	return value.value;
}

template <>
inline bool ECFValueHolder<std::string>::setValue(const std::string &value)
{
	this->value = value;
	return true;
}

template <>
inline bool ECFValueHolder<ECFExpression>::setValue(const std::string &value)
{
	this->value.value = value;
	return true;
}

template <>
inline bool ECFValueHolder<bool>::setValue(const std::string &value)
{
	if (value.size() == 0) {
		this->value = true;
		return true;
	} else {
		if (StringCompare::caseInsensitiveEq(value, "FALSE")) {
			this->value = false;
			return true;
		}
		if (StringCompare::caseInsensitiveEq(value, "TRUE")) {
			this->value = true;
			return true;
		}
		std::stringstream ss(value);
		ss >> this->value;
		return ss.eof();
	}
}

template <typename Ttype>
struct ECFEnumHolder: public ECFValue {
	Ttype &value;

	ECFEnumHolder(Ttype &value): value(value) {}

	std::string getValue() const { return metadata.options[static_cast<int>(value)].name; }

	bool setValue(const std::string &value)
	{
		for (size_t i = 0; i < metadata.options.size(); i++) {
			if (StringCompare::caseInsensitiveEq(value, metadata.options[i].name)) {
				this->value = static_cast<Ttype>(i);
				return true;
			}
		}
		size_t index;
		if (ECFValueHolder<size_t>(index).setValue(value)) {
			this->value = static_cast<Ttype>(index);
			return true;
		}
		return false;
	}

};

}



#endif /* SRC_CONFIG_VALUEHOLDER_H_ */
