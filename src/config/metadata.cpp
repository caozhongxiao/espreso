
#include "metadata.h"
#include "conditions.h"

#include "basis/logging/logging.h"

using namespace espreso;

std::string SIUnit::unit() const
{
	auto exponent = [] (const std::string &value, int exponent) {
		if (exponent == 0) {
			return std::string();
		}
		return value + std::to_string(exponent);
	};

	return
			exponent("m", metre) +
			exponent("kg", kilogram) +
			exponent("s", second) +
			exponent("A", ampere) +
			exponent("K", kelvin) +
			exponent("mol", mole) +
			exponent("cd", candela);
}

void ECFMetaData::checkdescription(const std::string &name, size_t size) const
{
	if (description.size() != size) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: " << name << " has incorrect number of descriptions.";
	}
}

void ECFMetaData::checkdatatype(const std::string &name, size_t size) const
{
	if (datatype.size() != size) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: " << name << " has incorrect number of datatypes.";
	}
}

void ECFMetaData::checkpattern(const std::string &name, size_t size) const
{
	if (pattern.size() != size) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: " << name << " has incorrect number of pattern values.";
	}
}

ECFMetaData& ECFMetaData::addconstraint(const ECFCondition &condition)
{
	delete this->condition; this->condition = condition.copy(); return *this;
}

template <typename TType>
static std::vector<TType> getsuffix(size_t start, const std::vector<TType> &data)
{
	if (start < data.size()) {
		return std::vector<TType>(data.begin() + start, data.end());
	}
	return std::vector<TType>();
}

ECFMetaData ECFMetaData::suffix(size_t start) const
{
	ECFMetaData ret(*this);
	ret.setdescription(getsuffix(start, description));
	ret.setdatatype(getsuffix(start, datatype));
	ret.setpattern(getsuffix(start, pattern));
	ret.tensor = NULL;
	ret.regionMap = NULL;
	return ret;
}

ECFMetaData::ECFMetaData(): tensor(NULL), regionMap(NULL)
{
	condition = new ECFCondition();
	isallowed = [] () { return true; };
	ismandatory = [] () { return true; };
}

ECFMetaData::ECFMetaData(const ECFMetaData &other)
: name(other.name), description(other.description),
  datatype(other.datatype), pattern(other.pattern),
  options(other.options), variables(other.variables),
  tensor(other.tensor), regionMap(other.regionMap),
  unit(other.unit), condition(other.condition->copy()),
  isallowed(other.isallowed), ismandatory(other.ismandatory)
{

}

ECFMetaData::ECFMetaData(ECFMetaData &&other)
: name(std::move(other.name)), description(std::move(other.description)),
  datatype(std::move(other.datatype)), pattern(std::move(other.pattern)),
  options(std::move(other.options)), variables(std::move(other.variables)),
  tensor(std::move(other.tensor)), regionMap(std::move(other.regionMap)),
  unit(std::move(other.unit)), condition(other.condition->copy()),
  isallowed(std::move(other.isallowed)), ismandatory(std::move(other.ismandatory))
{

}

ECFMetaData& ECFMetaData::operator=(const ECFMetaData &other)
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

ECFMetaData::~ECFMetaData()
{
	delete condition;
}

