
#ifndef SRC_CONFIG_OBJECTHOLDER_H_
#define SRC_CONFIG_OBJECTHOLDER_H_

#include "configuration.h"
#include "valueholder.h"

namespace espreso {

template <typename TParameter, typename TValue>
struct ECFValueMap: public ECFObject {
	std::map<TParameter, TValue> &value;

	ECFValueMap(std::map<TParameter, TValue> &value): value(value) {}

	virtual ECFParameter* getParameter(const std::string &name)
	{
		if (ECFObject::getParameter(name) != NULL) {
			return ECFObject::getParameter(name);
		}
		TParameter key;
		if (ECFValueHolder<TParameter>(key).setValue(name)) {
			return registerParameter(name, new ECFValueHolder<TValue>(value[key]), metadata.suffix(1));
		} else {
			return NULL;
		}
	}

	virtual void dropParameter(ECFParameter *parameter)
	{
		TParameter key;
		ECFValueHolder<TParameter>(key).setValue(parameter->name);
		value.erase(key);
		ECFObject::dropParameter(parameter);
	}

	virtual const ECFParameter* getPattern() const
	{
		ECFParameter *parameter = new ECFValueHolder<TValue>(_patternValue);
		parameter->setValue(metadata.pattern[1]);
		return registerPatternParameter(parameter);
	}

	virtual const void* data() const { return &value; }

private:
	mutable TValue _patternValue;
};

template <typename TParameter, typename TObject, typename... TArgs>
struct ECFObjectMap: public ECFObject {
	std::map<TParameter, TObject> &value;
	std::tuple<TArgs...> args;

	ECFObjectMap(std::map<TParameter, TObject> &value, TArgs... args): value(value), args(args...) {}
	ECFObjectMap(std::map<TParameter, TObject> &value, std::tuple<TArgs...> &args): value(value), args(args) {}

	virtual ECFParameter* getParameter(const std::string &name)
	{
		if (ECFObject::getParameter(name) != NULL) {
			return ECFObject::getParameter(name);
		}
		TParameter key;
		if (ECFValueHolder<TParameter>(key).setValue(name)) {
			auto it = value.emplace(std::piecewise_construct, std::forward_as_tuple(key), args);
			parameters.push_back(&it.first->second);
			parameters.back()->name = name;
			parameters.back()->metadata = metadata.suffix(1);
			return parameters.back();
		} else {
			return NULL;
		}
	}

	virtual void dropParameter(ECFParameter *parameter)
	{
		TParameter key;
		ECFValueHolder<TParameter>(key).setValue(parameter->name);
		value.erase(key);
		ECFObject::dropParameter(parameter);
	}

	virtual const ECFParameter* getPattern() const
	{
		std::map<TParameter, TObject*> dummy;
		auto it = value.emplace(std::piecewise_construct, std::forward_as_tuple(TParameter{}), args);
		return registerPatternParameter(&it.first->second);
	}

	virtual const void* data() const { return &value; }
};

template <typename TParameter1, typename TParameter2, typename TValue>
struct ECFValueMapMap: public ECFObject {
	std::map<TParameter1, std::map<TParameter2, TValue> > &value;

	ECFValueMapMap(std::map<TParameter1, std::map<TParameter2, TValue> > &value): value(value) {}

	virtual ECFParameter* getParameter(const std::string &name)
	{
		if (ECFObject::getParameter(name) != NULL) {
			return ECFObject::getParameter(name);
		}
		TParameter1 key;
		if (ECFValueHolder<TParameter1>(key).setValue(name)) {
			return registerParameter(name, new ECFValueMap<TParameter2, TValue>(value[key]), metadata.suffix(1));
		} else {
			return NULL;
		}
	}

	virtual void dropParameter(ECFParameter *parameter)
	{
		TParameter1 key;
		ECFValueHolder<TParameter1>(key).setValue(parameter->name);
		value.erase(key);
		ECFObject::dropParameter(parameter);
	}

	virtual const ECFParameter* getPattern() const
	{
		return registerPatternParameter(new ECFValueMap<TParameter2, TValue>(_patternValue));
	}

	virtual const void* data() const { return &value; }

private:
	mutable std::map<TParameter2, TValue> _patternValue;
};


template <typename TParameter1, typename TParameter2, typename TObject, typename... TArgs>
struct ECFObjectMapMap: public ECFObject {
	std::map<TParameter1, std::map<TParameter2, TObject> > &value;
	std::tuple<TArgs...> args;

	ECFObjectMapMap(std::map<TParameter1, std::map<TParameter2, TObject> > &value, TArgs... args): value(value), args(args...) {}

	virtual ECFParameter* getParameter(const std::string &name)
	{
		if (ECFObject::getParameter(name) != NULL) {
			return ECFObject::getParameter(name);
		}
		TParameter1 key;
		if (ECFValueHolder<TParameter1>(key).setValue(name)) {
			return registerParameter(name, new ECFObjectMap<TParameter2, TObject, TArgs...>(value[key], args), metadata.suffix(1));
		} else {
			return NULL;
		}
	}

	virtual void dropParameter(ECFParameter *parameter)
	{
		TParameter1 key;
		ECFValueHolder<TParameter1>(key).setValue(parameter->name);
		value.erase(key);
		ECFObject::dropParameter(parameter);
	}

	virtual const ECFParameter* getPattern() const
	{
		return registerPatternParameter(new ECFObjectMap<TParameter2, TObject>(_patternValue));
	}

	virtual const void* data() const { return &value; }

private:
	mutable std::map<TParameter2, TObject> _patternValue;
};

}

#endif /* SRC_CONFIG_OBJECTHOLDER_H_ */
