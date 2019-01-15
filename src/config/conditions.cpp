
#include "conditions.h"
#include "configuration.h"

using namespace espreso;


std::string ECFCondition::compose(const ECFParameter* parameter) const
{
	switch (parameter->metadata.datatype.front()) {
	case ECFDataType::ENUM_FLAGS:
	case ECFDataType::OPTION:
		return parameter->name + "=" + parameter->metadata.options[dynamic_cast<const EnumValue*>(value)->index()].name;
	default:
		return parameter->name + "=" + value->tostring();
	}
}

