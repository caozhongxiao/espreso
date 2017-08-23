
#ifndef SRC_CONFIGURATION_CONFIGURATION_HPP_
#define SRC_CONFIGURATION_CONFIGURATION_HPP_

#include <type_traits>
#include <iostream>

#include "configuration.h"
#include "valueholder.h"
#include "objectholder.h"

#define REGISTER(parameter, metadata) registerParameter(#parameter, parameter, metadata)
#define PNAME(parameter) #parameter

namespace espreso {

/////////// PARAMETER ///////////
/////////////////////////////////

// CLASS - STD::STRING (has to inherit from ECFObject)
template <typename Ttype>
typename std::enable_if<std::is_class<Ttype>::value && !std::is_same<Ttype, std::string>::value && !std::is_same<Ttype, ECFExpression>::value>::type
ECFObject::registerParameter(const std::string &name, Ttype &parameter, const ECFMetaData &metadata)
{
	metadata.checkdescription(name, 1);
	metadata.checkdatatype(name, 0);
	metadata.checkdefaultvalue(name, 0);

	parameters.push_back(&parameter);
	parameters.back()->name = name;
	parameters.back()->metadata = metadata;
}

// ENUM
template<typename Ttype>
typename std::enable_if<std::is_enum<Ttype>::value>::type
ECFObject::registerParameter(const std::string &name, Ttype &parameter, const ECFMetaData &metadata)
{
	metadata.checkdescription(name, 1);
	metadata.checkdatatype(name, 1);
	metadata.checkdefaultvalue(name, 0);

	registerParameter(name, new ECFEnumHolder<Ttype>(parameter), metadata);
}

// REST
template <typename Ttype>
typename std::enable_if<(!std::is_class<Ttype>::value && !std::is_enum<Ttype>::value) || std::is_same<Ttype, std::string>::value || std::is_same<Ttype, ECFExpression>::value>::type
ECFObject::registerParameter(const std::string &name, Ttype &parameter, const ECFMetaData &metadata)
{
	metadata.checkdescription(name, 1);
	metadata.checkdatatype(name, 1);
	metadata.checkdefaultvalue(name, 0);

	registerParameter(name, new ECFValueHolder<Ttype>(parameter), metadata);
}

////////////// MAP //////////////
/////////////////////////////////

// TYPE2 = CLASS - STD::STRING (has to inherit from ECFObject)
template<typename Ttype1, typename Ttype2>
typename std::enable_if<std::is_class<Ttype2>::value && !std::is_same<Ttype2, std::string>::value>::type
ECFObject::registerParameter(const std::string &name, std::map<Ttype1, Ttype2> &parameter, const ECFMetaData &metadata)
{
	metadata.checkdescription(name, 2);
	metadata.checkdatatype(name, 1);
	metadata.checkdefaultvalue(name, 1);

	registerParameter(name, new ECFObjectMap<Ttype1, Ttype2>(parameter), metadata);
}

// TYPE2 = REST
template<typename Ttype1, typename Ttype2>
typename std::enable_if<!std::is_class<Ttype2>::value || std::is_same<Ttype2, std::string>::value>::type
ECFObject::registerParameter(const std::string &name, std::map<Ttype1, Ttype2> &parameter, const ECFMetaData &metadata)
{
	metadata.checkdescription(name, 2);
	metadata.checkdatatype(name, 2);
	metadata.checkdefaultvalue(name, 2);

	registerParameter(name, new ECFValueMap<Ttype1, Ttype2>(parameter), metadata);
}

//////////// MAP MAP ////////////
/////////////////////////////////

// TYPE3 = CLASS - STD::STRING (has to inherit from ECFObject)
template<typename Ttype1, typename Ttype2, typename Ttype3>
typename std::enable_if<std::is_class<Ttype3>::value && !std::is_same<Ttype3, std::string>::value>::type
ECFObject::registerParameter(const std::string &name, std::map<Ttype1, std::map<Ttype2, Ttype3> > &parameter, const ECFMetaData &metadata)
{
	metadata.checkdescription(name, 3);
	metadata.checkdatatype(name, 2);
	metadata.checkdefaultvalue(name, 2);

	registerParameter(name, new ECFObjectMapMap<Ttype1, Ttype2, Ttype3>(parameter), metadata);
}

// TYPE3 = REST
template<typename Ttype1, typename Ttype2, typename Ttype3>
typename std::enable_if<!std::is_class<Ttype3>::value || std::is_same<Ttype3, std::string>::value>::type
ECFObject::registerParameter(const std::string &name, std::map<Ttype1, std::map<Ttype2, Ttype3> > &parameter, const ECFMetaData &metadata)
{
	metadata.checkdescription(name, 3);
	metadata.checkdatatype(name, 3);
	metadata.checkdefaultvalue(name, 3);

	registerParameter(name, new ECFValueMapMap<Ttype1, Ttype2, Ttype3>(parameter), metadata);
}

}


#endif /* SRC_CONFIGURATION_CONFIGURATION_HPP_ */
