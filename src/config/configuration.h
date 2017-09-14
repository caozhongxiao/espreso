
#ifndef SRC_CONFIG_CONFIGURATION_H_
#define SRC_CONFIG_CONFIGURATION_H_

#include <string>
#include <vector>
#include <map>
#include <functional>

namespace espreso {

struct Point;
class Evaluator;
struct TensorConfiguration;

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
	MATERIAL,
	LOAD_STEP,
	EXPRESSION,
	TENSOR,
	INTERVAL,
	SPACE,
	SEPARATOR
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

struct ECFExpression {
	std::string value;
	Evaluator *evaluator;

	double evaluate(const Point &p, double time = 0, double temperature = 0, double pressure = 0, double velocity = 0) const;
	ECFExpression();
	~ECFExpression();
};

struct ECFMetaData {
	std::vector<std::string> description;
	std::vector<ECFDataType> datatype;
	std::vector<std::string> pattern;
	std::vector<ECFOption> options;
	std::vector<std::string> variables;
	TensorConfiguration *tensor;
	std::string unit;

	std::function<bool(void)> isallowed;

	ECFMetaData& setdescription(const std::vector<std::string> &description) { this->description = description; return *this; }
	ECFMetaData& setdatatype(const std::vector<ECFDataType> &datatype) { this->datatype = datatype; return *this; }
	ECFMetaData& setpattern(const std::vector<std::string> &pattern) { this->pattern = pattern; return *this; }
	ECFMetaData& setvariables(const std::vector<std::string> &variables) { this->variables = variables; return *this; }
	ECFMetaData& settensor(TensorConfiguration &tensor) { this->tensor = &tensor; return *this; }
	ECFMetaData& setunit(const std::string &unit) { this->unit = unit; return *this; }
	ECFMetaData& allowonly(std::function<bool(void)> isallowed) { this->isallowed = isallowed; return *this; }

	ECFMetaData& addoption(const ECFOption &option) { options.push_back(option); return *this; }

	ECFMetaData& setcoordinatevariables() { return setvariables({ "X", "Y", "Z" }); }
	ECFMetaData& setboundaryconditionvariables() { return setvariables({ "X", "Y", "Z", "TIME" }); }
	ECFMetaData& setmaterialvariables() { return setvariables({ "X", "Y", "Z", "TIME", "TEMPERATURE" }); }

	void checkdescription(const std::string &name, size_t size) const;
	void checkdatatype(const std::string &name, size_t size) const;
	void checkpattern(const std::string &name, size_t size) const;

	ECFMetaData suffix(size_t start) const;

	ECFMetaData(): tensor(NULL) { isallowed = [] () { return true; }; }
};

struct ECFParameter {
	enum class Event {
		VALUE_SET
	};

	std::string name;
	ECFMetaData metadata;

	virtual bool isValue() const =0;
	virtual bool isObject() const =0;

	virtual bool setValue(const std::string &value);
	virtual std::string getValue() const =0;
	virtual ECFParameter* getParameter(const std::string &name) =0;
	virtual ECFParameter* getParameter(const void* data) =0;

	virtual const ECFParameter* getPattern() const =0;
	virtual const void* data() const =0;

	virtual void addListener(Event event, std::function<void()> listener);

	virtual ~ECFParameter() {};
protected:
	virtual bool _setValue(const std::string &value) =0;
	std::vector<std::function<void()> > _setValueListeners;
};

struct ECFSeparator: public ECFParameter {

	bool isValue() const { return false; }
	bool isObject() const { return false; }

	bool _setValue(const std::string &value) { return false; }
	std::string getValue() const { return ""; }
	ECFParameter* getParameter(const std::string &name) { return NULL; }
	ECFParameter* getParameter(const void* data) { return NULL; }

	virtual const ECFParameter* getPattern() const { return NULL; }
	virtual const void* data() const { return NULL; }
};

struct ECFValue: public ECFParameter {

	virtual ECFParameter* getParameter(const std::string &name) { return NULL; }
	virtual ECFParameter* getParameter(const void* data) { return NULL; }

	bool isValue() const { return true; }
	bool isObject() const { return false; }

	virtual const ECFParameter* getPattern() const { return NULL; }
};

struct ECFObject: public ECFParameter {
	std::vector<ECFParameter*> parameters;

	bool isValue() const { return false; }
	bool isObject() const { return true; }

	virtual bool _setValue(const std::string &value);
	virtual std::string getValue() const;
	virtual ECFParameter* getParameter(const std::string &name);
	virtual ECFParameter* getParameter(const void* data);

	virtual const ECFParameter* getPattern() const { return NULL; }
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


	virtual ~ECFObject();

	void dropParameter(ECFParameter *parameter);
protected:
	ECFParameter* addSeparator();
	ECFParameter* addSpace();

	/////////// PARAMETER ///////////
	/////////////////////////////////

	// CLASS - STD::STRING - ECFExpression (has to inherit from ECFObject)
	template<typename Ttype>
	typename std::enable_if<std::is_class<Ttype>::value && !std::is_same<Ttype, std::string>::value && !std::is_same<Ttype, ECFExpression>::value, ECFParameter*>::type
	registerParameter(const std::string &name, Ttype &parameter, const ECFMetaData &metadata);

	// ENUM
	template<typename Ttype>
	typename std::enable_if<std::is_enum<Ttype>::value, ECFParameter*>::type
	registerParameter(const std::string &name, Ttype &parameter, const ECFMetaData &metadata);

	// REST
	template<typename Ttype>
	typename std::enable_if<(!std::is_class<Ttype>::value && !std::is_enum<Ttype>::value) || std::is_same<Ttype, std::string>::value || std::is_same<Ttype, ECFExpression>::value, ECFParameter*>::type
	registerParameter(const std::string &name, Ttype &parameter, const ECFMetaData &metadata);


	////////////// MAP //////////////
	/////////////////////////////////

	// TYPE2 = CLASS - STD::STRING (has to inherit from ECFObject)
	template<typename Ttype1, typename Ttype2, typename... TArgs>
	typename std::enable_if<std::is_class<Ttype2>::value && !std::is_same<Ttype2, std::string>::value, ECFParameter*>::type
	registerParameter(const std::string &name, std::map<Ttype1, Ttype2> &parameter, const ECFMetaData &metadata, TArgs... args);

	// TYPE2 = REST
	template<typename Ttype1, typename Ttype2>
	typename std::enable_if<!std::is_class<Ttype2>::value || std::is_same<Ttype2, std::string>::value, ECFParameter*>::type
	registerParameter(const std::string &name, std::map<Ttype1, Ttype2> &parameter, const ECFMetaData &metadata);

	//////////// MAP MAP ////////////
	/////////////////////////////////

	// TYPE3 = CLASS - STD::STRING (has to inherit from ECFObject)
	template<typename Ttype1, typename Ttype2, typename Ttype3, typename... TArgs>
	typename std::enable_if<std::is_class<Ttype3>::value && !std::is_same<Ttype3, std::string>::value, ECFParameter*>::type
	registerParameter(const std::string &name, std::map<Ttype1, std::map<Ttype2, Ttype3> > &parameter, const ECFMetaData &metadata, TArgs... args);

	// TYPE3 = REST
	template<typename Ttype1, typename Ttype2, typename Ttype3>
	typename std::enable_if<!std::is_class<Ttype3>::value || std::is_same<Ttype3, std::string>::value, ECFParameter*>::type
	registerParameter(const std::string &name, std::map<Ttype1, std::map<Ttype2, Ttype3> > &parameter, const ECFMetaData &metadata);


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
