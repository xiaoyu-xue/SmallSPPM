#pragma once
#include "common/Core.h"

GY_NAMESPACE_BEGIN

template<typename T>
T* CreateInstanceRawPtr(std::string alias);

template<typename T>
std::shared_ptr<T> CreateInstanceSharedPtr(std::string alias);

template<typename T>
std::unique_ptr<T> CreateInstanceUniquePtr(std::string alias);


template<typename T>
class FactoryBase
{
	using CreateMethod = std::function<T* ()>;
protected:
	std::map<std::string, CreateMethod> mCreateMethods;
public:
	T* CreateInstanceRawPtr(std::string alias)
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

#define GY_INTERFACE(BaseClassName)															\
class GY_FACTORY_NAME(BaseClassName) : public FactoryBase<BaseClassName>					\
{																							\
																							\
};																							


#define GY_INTERFACE_DEF(BaseClassName)														\
template<>																					\
BaseClassName* CreateInstanceRawPtr<BaseClassName>(std::string alias)						\
{																							\
	auto factory = GY_FACTORY_NAME(BaseClassName)::GetInstancePtr();						\
	return factory->CreateInstanceRawPtr(alias);											\
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