
#ifndef SRC_CONFIG_CONDITIONS_H_
#define SRC_CONFIG_CONDITIONS_H_

#include <string>

namespace espreso {

class ECFParameter;

struct GeneralValue {
	virtual bool equal(const void *data) const { return true; }
	virtual std::string tostring() const { return "TRUE"; }

	virtual ~GeneralValue() {}
	virtual GeneralValue* copy() const { return new GeneralValue(); }
};

struct EnumValue: public GeneralValue {
	virtual int index() const =0;
};

template <typename TValue>
struct EnumValueHolder: public EnumValue {
	TValue value;

	bool equal(const void *data) const { return *static_cast<const TValue*>(data) == value; }
	int index() const { return static_cast<int>(value); }

	EnumValueHolder(const TValue &value): value(value) {}
	GeneralValue* copy() const { return new EnumValueHolder<TValue>(value); }
};

template <typename TValue>
struct GeneralValueHolder: public GeneralValue {
	TValue value;

	bool equal(const void *data) const { return *static_cast<const TValue*>(data) == value; }
	std::string tostring() const { return std::to_string(value); }

	GeneralValueHolder(const TValue &value): value(value) {}
	GeneralValue* copy() const { return new GeneralValueHolder<TValue>(value); }
};

class ECFCondition {

protected:
	const void *parameter;
	const GeneralValue *value;

public:
	virtual bool evaluate() const { return !parameter || value->equal(parameter); }

	bool isset() const { return parameter; }
	virtual bool match(const void* parameter) const { return this->parameter == parameter; }
	virtual std::string compose(const ECFParameter* parameter) const;

	ECFCondition()
	: parameter(NULL), value(new GeneralValue()) {}

	template <typename TValue>
	typename std::enable_if<std::is_enum<TValue>::value, const GeneralValue*>::type
	init(const TValue &value) { return new EnumValueHolder<TValue>(value); }

	template <typename TValue>
	typename std::enable_if<!std::is_enum<TValue>::value, const GeneralValue*>::type
	init(const TValue &value) { return new GeneralValueHolder<TValue>(value); }

	template <typename TParemeter, typename TValue>
	ECFCondition(const TParemeter &parameter, const TValue &value)
	: parameter(&parameter), value(init<TValue>(value)) {}

	ECFCondition(const ECFCondition &other)
	: parameter(other.parameter), value(other.value->copy()) {}

	ECFCondition(ECFCondition &&other)
	: parameter(std::move(other.parameter)), value(other.value->copy()) {}

	virtual ECFCondition* copy() const { return new ECFCondition(*this); }

	virtual ~ECFCondition() { delete value; }
};

class ECFNotCondition: public ECFCondition {

public:
	virtual bool evaluate() const { return !parameter || !value->equal(parameter); }

	ECFNotCondition() {}

	template <typename TParemeter, typename TValue>
	ECFNotCondition(const TParemeter &parameter, const TValue &value)
	: ECFCondition(parameter, value) {}

	ECFNotCondition(const ECFCondition &other)
	: ECFCondition(other) {}

	ECFNotCondition(const ECFNotCondition &other)
	: ECFCondition(other) {}

	ECFNotCondition(ECFCondition &&other)
	: ECFCondition(std::move(other)) {}

	virtual ECFCondition* copy() const { return new ECFNotCondition(*this); }
};

inline ECFNotCondition operator!(const ECFCondition &other) { return other; }

}

#endif /* SRC_CONFIG_CONDITIONS_H_ */
