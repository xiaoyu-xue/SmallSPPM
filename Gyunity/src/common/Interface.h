#pragma once
#include "common/Core.h"

GY_NAMESPACE_BEGIN

template<typename T>
inline T* CreateInstanceRawPtr(std::string alias);


class InstanceCreator
{
public:
	template<typename T>
	inline T* CreateRawPtr(std::string alias)
	{
		return nullptr;
	}
	static InstanceCreator* GetInstancePtr()
	{
		static InstanceCreator factoryHolder;
		return &factoryHolder;
	}
};

template<typename T>
class FactoryBase
{
	using CreateMethod = std::function<T* ()>;
protected:
	std::map<std::string, CreateMethod> mCreateMethods;
public:
	T* CreateRawPtr(std::string alias)
	{
		if (mCreateMethods.find(alias) != mCreateMethods.end())
			return mCreateMethods[alias]();
		return nullptr;
	}
	void Register(std::string alias, CreateMethod createMethod)
	{
		mCreateMethods[alias] = createMethod;
	}
	static FactoryBase<T>* GetInstancePtr()
	{
		static FactoryBase<T> instance;
		return &instance;
	}
};


#define GY_FACTORY_NAME(T) Factory_##T
#define GY_IMPLEMENTATION_NAME(T) Impementation_##T##_Injector

#define GY_INTERFACE_DEF(BaseClassName)														\
class GY_FACTORY_NAME(BaseClassName) : public FactoryBase<BaseClassName>					\
{																							\
																							\
};																							\
template<>																					\
inline BaseClassName* CreateInstanceRawPtr<BaseClassName>(std::string alias)				\
{																							\
	auto factory = GY_FACTORY_NAME(BaseClassName)::GetInstancePtr();						\
	return factory->CreateRawPtr(alias);													\
}

#define GY_IMPLEMENTATION_DEF(BaseClassName, ClassName, Alias)								\
class GY_IMPLEMENTATION_NAME(ClassName)														\
{																							\
	using CreateMethod = std::function<BaseClassName*()>;									\
public:																						\
	GY_IMPLEMENTATION_NAME(ClassName)()														\
	{																						\
		auto factory = GY_FACTORY_NAME(BaseClassName)::GetInstancePtr();					\
		CreateMethod createMethod = []()->BaseClassName* {return new ClassName(); };		\
		factory->Register(Alias, createMethod);												\
	}																						\
}Impementation_##BaseClassName##_Instance_##ClassName;


GY_NAMESPACE_END