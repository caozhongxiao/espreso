
#ifndef SRC_CONFIGURATION_CONFIGURATION_HPP_
#define SRC_CONFIGURATION_CONFIGURATION_HPP_

#include <type_traits>
#include <iostream>

#include "configuration.h"
#include "valueholder.h"
#include "objectholder.h"

#define REGISTER(parameter, metadata, ...) registerParameter(#parameter, parameter, metadata, ##__VA_ARGS__)
#define PNAME(parameter) #parameter

namespace espreso {

/////////// PARAMETER ///////////
/////////////////////////////////

// CLASS - STD::STRING - ECFExpression (has to inherit from ECFObject)
template <typename Ttype>
typename std::enable_if<std::is_class<Ttype>::value && !std::is_same<Ttype, std::string>::value && !std::is_same<Ttype, ECFExpression>::value, ECFParameter*>::type
ECFObject::registerParameter(const std::string &name, Ttype &parameter, const ECFMetaData &metadata)
{
	metadata.checkdescription(name, 1);
	metadata.checkdatatype(name, 0);
	metadata.checkpattern(name, 0);

	parameters.push_back(&parameter);
	parameters.back()->name = name;
	parameters.back()->metadata = metadata;
	return parameters.back();
}

// ENUM
template<typename Ttype>
typename std::enable_if<std::is_enum<Ttype>::value, ECFParameter*>::type
ECFObject::registerParameter(const std::string &name, Ttype &parameter, const ECFMetaData &metadata)
{
	metadata.checkdescription(name, 1);
	metadata.checkdatatype(name, 1);
	metadata.checkpattern(name, 0);

	return registerParameter(name, new ECFEnumHolder<Ttype>(parameter), metadata);
}

// REST
template <typename Ttype>
typename std::enable_if<(!std::is_class<Ttype>::value && !std::is_enum<Ttype>::value) || std::is_same<Ttype, std::string>::value || std::is_same<Ttype, ECFExpression>::value, ECFParameter*>::type
ECFObject::registerParameter(const std::string &name, Ttype &parameter, const ECFMetaData &metadata)
{
	metadata.checkdescription(name, 1);
	metadata.checkdatatype(name, 1);
	metadata.checkpattern(name, 0);

	return registerParameter(name, new ECFValueHolder<Ttype>(parameter), metadata);
}

////////////// MAP //////////////
/////////////////////////////////

// TYPE2 = CLASS - STD::STRING (has to inherit from ECFObject)
template<typename Ttype1, typename Ttype2, typename... TArgs>
typename std::enable_if<std::is_class<Ttype2>::value && !std::is_same<Ttype2, std::string>::value, ECFParameter*>::type
ECFObject::registerParameter(const std::string &name, std::map<Ttype1, Ttype2> &parameter, const ECFMetaData &metadata, TArgs... args)
{
	metadata.checkdescription(name, 2);
	metadata.checkdatatype(name, 1);
	metadata.checkpattern(name, 1);

	return registerParameter(name, new ECFObjectMap<Ttype1, Ttype2, TArgs...>(parameter, args...), metadata);
}

// TYPE2 = REST
template<typename Ttype1, typename Ttype2>
typename std::enable_if<!std::is_class<Ttype2>::value || std::is_same<Ttype2, std::string>::value, ECFParameter*>::type
ECFObject::registerParameter(const std::string &name, std::map<Ttype1, Ttype2> &parameter, const ECFMetaData &metadata)
{
	metadata.checkdescription(name, 2);
	metadata.checkdatatype(name, 2);
	metadata.checkpattern(name, 2);

	return registerParameter(name, new ECFValueMap<Ttype1, Ttype2>(parameter), metadata);
}

//////////// MAP MAP ////////////
/////////////////////////////////

// TYPE3 = CLASS - STD::STRING (has to inherit from ECFObject)
template<typename Ttype1, typename Ttype2, typename Ttype3, typename... TArgs>
typename std::enable_if<std::is_class<Ttype3>::value && !std::is_same<Ttype3, std::string>::value, ECFParameter*>::type
ECFObject::registerParameter(const std::string &name, std::map<Ttype1, std::map<Ttype2, Ttype3> > &parameter, const ECFMetaData &metadata, TArgs... args)
{
	metadata.checkdescription(name, 3);
	metadata.checkdatatype(name, 2);
	metadata.checkpattern(name, 2);

	return registerParameter(name, new ECFObjectMapMap<Ttype1, Ttype2, Ttype3, TArgs...>(parameter, args...), metadata);
}

// TYPE3 = REST
template<typename Ttype1, typename Ttype2, typename Ttype3>
typename std::enable_if<!std::is_class<Ttype3>::value || std::is_same<Ttype3, std::string>::value, ECFParameter*>::type
ECFObject::registerParameter(const std::string &name, std::map<Ttype1, std::map<Ttype2, Ttype3> > &parameter, const ECFMetaData &metadata)
{
	metadata.checkdescription(name, 3);
	metadata.checkdatatype(name, 3);
	metadata.checkpattern(name, 3);

	return registerParameter(name, new ECFValueMapMap<Ttype1, Ttype2, Ttype3>(parameter), metadata);
}

}


#endif /* SRC_CONFIGURATION_CONFIGURATION_HPP_ */
