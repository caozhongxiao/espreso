
#ifndef SRC_CONFIG_CONFIGURATION_H_
#define SRC_CONFIG_CONFIGURATION_H_

#include <string>
#include <vector>
#include <map>
#include <functional>

#include "regionmap.h"

namespace espreso {

struct Point;
struct TensorConfiguration;
struct ECFParameter;

enum class ECFDataType {
	BOOL,
	STRING,
	INTEGER,
	POSITIVE_INTEGER,
	NONNEGATIVE_INTEGER,
	FLOAT,
	ENUM_FLAGS,
	OPTION,
	REGION,
	BOUNDARY_REGION,
	ELEMENTS_REGION,
	MATERIAL,
	LOAD_STEP,
	EXPRESSION,
	TENSOR,
	INTERVAL,
	SPACE,
	SEPARATOR
};

enum class DIMENSION {
	D1,
	D2,
	D3,
	Z
};

struct ECFOption {
	std::string name;
	std::string description;
	std::function<bool(void)> isallowed;

	ECFOption& setname(const std::string &name) { this->name = name; return *this; }
	ECFOption& setdescription(const std::string &description) { this->description = description; return *this; }
	ECFOption& allowonly(std::function<bool(void)> isallowed) { this->isallowed = isallowed; return *this; }

	ECFOption() { isallowed = [] () { return true; }; }
};

struct SIUnit {
	int metre, kilogram, second, ampere, kelvin, mole, candela;

	SIUnit(): metre(0), kilogram(0), second(0), ampere(0), kelvin(0), mole(0), candela(0) {}

	SIUnit(int metre, int kilogram, int second, int ampere, int kelvin, int mole, int candela)
	: metre(metre), kilogram(kilogram), second(second), ampere(ampere), kelvin(kelvin), mole(mole), candela(candela)
	{}

	std::string unit() const;
};

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

struct ECFMetaData {
	std::string name;
	std::vector<std::string> description;
	std::vector<ECFDataType> datatype;
	std::vector<std::string> pattern;
	std::vector<ECFOption> options;
	std::vector<std::string> variables;
	TensorConfiguration *tensor;
	RegionMapBase *regionMap;
	SIUnit unit;

	ECFCondition *condition;

	std::function<bool(void)> isallowed;
	std::function<bool(void)> ismandatory;

	ECFMetaData& setname(const std::string &name) { this->name = name; return *this; }
	ECFMetaData& setdescription(const std::vector<std::string> &description) { this->description = description; return *this; }
	ECFMetaData& setdatatype(const std::vector<ECFDataType> &datatype) { this->datatype = datatype; return *this; }
	ECFMetaData& setpattern(const std::vector<std::string> &pattern) { this->pattern = pattern; return *this; }
	ECFMetaData& setvariables(const std::vector<std::string> &variables) { this->variables = variables; return *this; }
	ECFMetaData& settensor(TensorConfiguration &tensor) { this->tensor = &tensor; return *this; }
	ECFMetaData& setRegionMap(RegionMapBase &rMap) { this->regionMap = &rMap; return *this; }
	ECFMetaData& setunit(const SIUnit &unit) { this->unit = unit; return *this; }
	ECFMetaData& allowonly(std::function<bool(void)> isallowed) { this->isallowed = isallowed; return *this; }
	ECFMetaData& mandatoryonly(std::function<bool(void)> ismandatory) { this->ismandatory = ismandatory; return *this; }
	ECFMetaData& addconstraint(const ECFCondition &condition) { delete this->condition; this->condition = condition.copy(); return *this; }

	ECFMetaData& addoption(const ECFOption &option) { options.push_back(option); return *this; }

	static std::vector<std::string> getcoordinatevariables() { return { "X", "Y", "Z" }; }
	static std::vector<std::string> getboundaryconditionvariables() { return { "X", "Y", "Z", "TIME" }; }
	static std::vector<std::string> getmaterialvariables() { return { "X", "Y", "Z", "TIME", "TEMPERATURE" }; }

	void checkdescription(const std::string &name, size_t size) const;
	void checkdatatype(const std::string &name, size_t size) const;
	void checkpattern(const std::string &name, size_t size) const;

	ECFMetaData suffix(size_t start) const;

	ECFMetaData(): tensor(NULL), regionMap(NULL)
	{
		condition = new ECFCondition();
		isallowed = [] () { return true; };
		ismandatory = [] () { return true; };
	}

	ECFMetaData(const ECFMetaData &other)
	: name(other.name), description(other.description),
	  datatype(other.datatype), pattern(other.pattern),
	  options(other.options), variables(other.variables),
	  tensor(other.tensor), regionMap(other.regionMap),
	  unit(other.unit), isallowed(other.isallowed),
	  ismandatory(other.ismandatory), condition(other.condition->copy())
	{

	}

	ECFMetaData(ECFMetaData &&other)
	: name(std::move(other.name)), description(std::move(other.description)),
	  datatype(std::move(other.datatype)), pattern(std::move(other.pattern)),
	  options(std::move(other.options)), variables(std::move(other.variables)),
	  tensor(std::move(other.tensor)), regionMap(std::move(other.regionMap)),
	  unit(std::move(other.unit)), isallowed(std::move(other.isallowed)),
	  ismandatory(std::move(other.ismandatory)), condition(other.condition->copy())
	{

	}

	ECFMetaData& operator=(const ECFMetaData &other)
	{
		if (this != &other) {
			name = other.name;
			description = other.description;
			datatype = other.datatype;
			pattern = other.pattern;
			options = other.options;
			variables = other.variables;
			tensor = other.tensor;
			regionMap = other.regionMap;
			unit = other.unit;
			isallowed = other.isallowed;
			ismandatory = other.ismandatory;
			condition = other.condition->copy();
		}
		return *this;
	}

	~ECFMetaData()
	{
		delete condition;
	}
};

struct ECFParameter {
	enum class Event {
		VALUE_SET,
		PARAMETER_GET
	};

	std::string name;
	ECFMetaData metadata;

	virtual bool isValue() const =0;
	virtual bool isObject() const =0;

	virtual bool setValue(const std::string &value);
	virtual std::string getValue() const =0;
	virtual ECFParameter* getParameter(const std::string &name);
	virtual ECFParameter* getParameter(const char* name);
	virtual ECFParameter* getParameter(const void* data);

	virtual ECFParameter* getPattern() const =0;
	virtual const void* data() const =0;

	virtual void addListener(Event event, std::function<void(const std::string &value)> listener);

	virtual void defaultName();
	virtual bool isvisible() { return metadata.condition->evaluate(); }
	virtual ECFParameter* registerAdditionalParameter(ECFParameter* parameter);

	virtual ~ECFParameter() {};
protected:
	virtual bool _setValue(const std::string &value) =0;
	virtual ECFParameter* _getParameter(const std::string &name) =0;
	virtual ECFParameter* _getParameter(const void* data) =0;
	virtual ECFParameter* _triggerParameterGet(ECFParameter* parameter);

	std::vector<std::function<void(const std::string &value)> > _setValueListeners;
	std::vector<std::function<void(const std::string &name)> > _parameterGetListeners;
};

struct ECFSeparator: public ECFParameter {

	bool isValue() const { return false; }
	bool isObject() const { return false; }

	std::string getValue() const { return ""; }

	virtual ECFParameter* getPattern() const { return NULL; }
	virtual const void* data() const { return NULL; }

protected:
	bool _setValue(const std::string &value) { return false; }
	ECFParameter* _getParameter(const std::string &name) { return NULL; }
	ECFParameter* _getParameter(const void* data) { return NULL; }
};

struct ECFValue: public ECFParameter {

	bool isValue() const { return true; }
	bool isObject() const { return false; }

	virtual ECFParameter* getPattern() const { return NULL; }

protected:
	ECFParameter* _getParameter(const std::string &name) { return NULL; }
	ECFParameter* _getParameter(const void* data) { return NULL; }
};

struct ECFObject: public ECFParameter {
	std::vector<ECFParameter*> parameters;

	bool isValue() const { return false; }
	bool isObject() const { return true; }

	virtual std::string getValue() const;

	virtual ECFParameter* getPattern() const { return NULL; }
	virtual const void* data() const { return this; }

	void forEachParameters(std::function<void(ECFParameter*)> fnc, bool onlyAllowed = true);
	void forEachParameters(std::function<void(const ECFParameter*)> fnc, bool onlyAllowed = true) const;

	ECFObject() {}

	// Assigning of parameters invalidates set/get methods -> skip it
	ECFObject& operator=(ECFObject &other) { return *this; }
	ECFObject& operator=(const ECFObject &other) { return *this; }

	// Copy constructor skip default constructor -> disable it
	ECFObject(ECFObject &other) = delete;
	ECFObject(const ECFObject &other) = delete;

	// Never use move constructors
	ECFObject& operator=(ECFObject &&other) = delete;
	ECFObject& operator=(const ECFObject &&other) = delete;
	ECFObject(ECFObject &&other) = delete;
	ECFObject(const ECFObject &&other) = delete;

	virtual ~ECFObject();

	virtual ECFParameter* registerAdditionalParameter(ECFParameter* parameter);
	virtual void dropParameter(ECFParameter *parameter);
	virtual void dropAllParameters();
protected:
	ECFParameter* addSeparator();
	ECFParameter* addSpace();

	virtual bool _setValue(const std::string &value);
	virtual ECFParameter* _getParameter(const std::string &name);
	virtual ECFParameter* _getParameter(const void* data);

	/////////// PARAMETER ///////////
	/////////////////////////////////

	// Child of ECFObject
	template<typename Ttype>
	typename std::enable_if<std::is_class<Ttype>::value && std::is_base_of<ECFObject, Ttype>::value, ECFParameter*>::type
	registerParameter(const std::string &name, Ttype &parameter, const ECFMetaData &metadata);

	// ENUM
	template<typename Ttype>
	typename std::enable_if<std::is_enum<Ttype>::value, ECFParameter*>::type
	registerParameter(const std::string &name, Ttype &parameter, const ECFMetaData &metadata);

	// REST
	template<typename Ttype>
	typename std::enable_if<(!std::is_class<Ttype>::value && !std::is_enum<Ttype>::value) || (std::is_class<Ttype>::value && !std::is_base_of<ECFObject, Ttype>::value), ECFParameter*>::type
	registerParameter(const std::string &name, Ttype &parameter, const ECFMetaData &metadata);


	////////////// MAP //////////////
	/////////////////////////////////

	// TYPE2 = Child of ECFObject
	template<typename Ttype1, typename Ttype2, typename... TArgs>
	typename std::enable_if<std::is_class<Ttype2>::value && std::is_base_of<ECFObject, Ttype2>::value, ECFParameter*>::type
	registerParameter(const std::string &name, std::map<Ttype1, Ttype2> &parameter, const ECFMetaData &metadata, TArgs... args);

	// TYPE2 = ENUM
	template<typename Ttype1, typename Ttype2>
	typename std::enable_if<std::is_enum<Ttype2>::value, ECFParameter*>::type
	registerParameter(const std::string &name, std::map<Ttype1, Ttype2> &parameter, const ECFMetaData &metadata);

	// TYPE2 = REST
	template<typename Ttype1, typename Ttype2, typename... TArgs>
	typename std::enable_if<(!std::is_class<Ttype2>::value && !std::is_enum<Ttype2>::value) || (std::is_class<Ttype2>::value && !std::is_base_of<ECFObject, Ttype2>::value), ECFParameter*>::type
	registerParameter(const std::string &name, std::map<Ttype1, Ttype2> &parameter, const ECFMetaData &metadata, TArgs... args);

	//////////// MAP MAP ////////////
	/////////////////////////////////

	// TYPE3 = Child of ECFObject
	template<typename Ttype1, typename Ttype2, typename Ttype3, typename... TArgs>
	typename std::enable_if<std::is_class<Ttype3>::value && std::is_base_of<ECFObject, Ttype3>::value, ECFParameter*>::type
	registerParameter(const std::string &name, std::map<Ttype1, std::map<Ttype2, Ttype3> > &parameter, const ECFMetaData &metadata, TArgs... args);

	// TYPE3 = REST
	template<typename Ttype1, typename Ttype2, typename Ttype3>
	typename std::enable_if<!std::is_class<Ttype3>::value || (std::is_class<Ttype3>::value && !std::is_base_of<ECFObject, Ttype3>::value), ECFParameter*>::type
	registerParameter(const std::string &name, std::map<Ttype1, std::map<Ttype2, Ttype3> > &parameter, const ECFMetaData &metadata);


	////////// REGION MAP ///////////
	/////////////////////////////////
	template<typename Ttype, typename... TArgs>
	ECFParameter* registerParameter(const std::string &name, RegionMap<Ttype> &parameter, ECFMetaData &metadata, TArgs... args);

	/////////////////////////////////
	ECFParameter* getWithError(const std::string &name);
	void moveLastBefore(const std::string &name);
	ECFParameter* registerParameter(const std::string &name, ECFParameter *parameter, const ECFMetaData &metadata);
	ECFParameter* registerPatternParameter(ECFParameter *parameter) const;

private:
	mutable std::vector<ECFParameter*> registeredParameters;
};

}



#endif /* SRC_CONFIG_CONFIGURATION_H_ */
