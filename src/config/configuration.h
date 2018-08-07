
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

struct ECFMetaData {
	std::vector<std::string> description;
	std::vector<ECFDataType> datatype;
	std::vector<std::string> pattern;
	std::vector<ECFOption> options;
	std::vector<std::string> variables;
	TensorConfiguration *tensor;
	RegionMapBase *regionMap;
	std::string unit;

	std::function<bool(void)> isallowed;
	std::function<bool(void)> ismandatory;

	ECFMetaData& setdescription(const std::vector<std::string> &description) { this->description = description; return *this; }
	ECFMetaData& setdatatype(const std::vector<ECFDataType> &datatype) { this->datatype = datatype; return *this; }
	ECFMetaData& setpattern(const std::vector<std::string> &pattern) { this->pattern = pattern; return *this; }
	ECFMetaData& setvariables(const std::vector<std::string> &variables) { this->variables = variables; return *this; }
	ECFMetaData& settensor(TensorConfiguration &tensor) { this->tensor = &tensor; return *this; }
	ECFMetaData& setRegionMap(RegionMapBase &rMap) { this->regionMap = &rMap; return *this; }
	ECFMetaData& setunit(const std::string &unit) { this->unit = unit; return *this; }
	ECFMetaData& allowonly(std::function<bool(void)> isallowed) { this->isallowed = isallowed; return *this; }
	ECFMetaData& mandatoryonly(std::function<bool(void)> ismandatory) { this->ismandatory = ismandatory; return *this; }

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
		isallowed = [] () { return true; };
		ismandatory = [] () { return true; };
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
