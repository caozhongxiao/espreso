
#ifndef SRC_BASIS_PARAMETERS_ENVELOPES_ENVELOPE_H_
#define SRC_BASIS_PARAMETERS_ENVELOPES_ENVELOPE_H_

#include "../parser.h"

namespace espreso {

template <typename TValue>
struct Option {
	std::string name;
	TValue value;
	std::string description;
};

struct Envelope {
	virtual bool set(const std::string &value) =0;
	virtual std::string get() const =0;

	virtual Envelope* copy() =0;
	virtual size_t options() const
	{
		return 0;
	}

	virtual std::string optionName(size_t index) const
	{
		return "";
	}

	virtual std::string optionDesc(size_t index) const
	{
		return "";
	}

	virtual ~Envelope() {};
};

}

#endif /* SRC_BASIS_PARAMETERS_ENVELOPES_ENVELOPE_H_ */