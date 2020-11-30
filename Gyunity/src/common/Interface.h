#pragma once
#include "common/Core.h"

GY_NAMESPACE_BEGIN

template<typename T>
T* CreateRawPtr(const std::string& alias);

template<typename T>
std::shared_ptr<T> CreateSharedPtr(const std::string& alias);

template<typename T>
std::unique_ptr<T> CreateUniquePtr(const std::string& alias);

template<typename T>
T* CreatePlacementPtr(const std::string& alias, void* placement);


class Object 
{
public:
	virtual ~Object() {}
};


template<typename T>
class FactoryBase
{
	using CreateRawPtrMethod = std::function<T* ()>;
	using CreateSharedPtrMethod = std::function<std::shared_ptr<T> ()>;
	using CreateUniquePtrMethod = std::function<std::unique_ptr<T> ()>;
	using CreatePlacementMethod = std::function<T* (void*)>;
protected:
	std::map<std::string, CreateRawPtrMethod> mCreateRawPtrMethods;
	std::map<std::string, CreateSharedPtrMethod> mCreateSharedPtrMethods;
	std::map<std::string, CreateUniquePtrMethod> mCreateUniquePtrMethods;
	std::map<std::string, CreatePlacementMethod> mCreatePlacementMethods;
public:
	T* CreateRawPtr(const std::string& alias)
	{
		auto createMethod = mCreateRawPtrMethods.find(alias);
		if (createMethod != mCreateRawPtrMethods.end())
			return (createMethod->second)();
		return nullptr;
	}

	std::shared_ptr<T> CreateSharedPtr(const std::string& alias)
	{
		auto createMethod = mCreateSharedPtrMethods.find(alias);
		if (createMethod != mCreateSharedPtrMethods.end())
			return (createMethod->second)();
		return nullptr;
	}

	std::unique_ptr<T> CreateUniquePtr(const std::string& alias)
	{
		auto createMethod = mCreateUniquePtrMethods.find(alias);
		if (createMethod != mCreateUniquePtrMethods.end())
			return (createMethod->second)();
		return nullptr;
	}

	T* CreatePlacementPtr(const std::string& alias, void* placement)
	{
		auto createMethod = mCreatePlacementMethods.find(alias);
		if (createMethod != mCreatePlacementMethods.end())
			return (createMethod->second)(placement);
		return nullptr;
	}

	template<typename ProductType>
	void Register(std::string alias)
	{
		mCreateRawPtrMethods.insert(
			std::make_pair(alias, []()->T* { return new ProductType(); })
		);
		mCreateSharedPtrMethods.insert(
			std::make_pair(alias, []()->std::shared_ptr<T> { return std::make_shared<ProductType>(); })
		);
		mCreateUniquePtrMethods.insert(
			std::make_pair(alias, []()->std::unique_ptr<T> { return std::make_unique<ProductType>(); })
		);
		mCreatePlacementMethods.insert(
			std::make_pair(alias, [](void* placement)->T* { return new (placement) ProductType(); })
		);

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
																							\
template<>																					\
inline BaseClassName* CreateRawPtr<BaseClassName>(const std::string &alias)					\
{																							\
	auto factory = GY_FACTORY_NAME(BaseClassName)::GetInstancePtr();						\
	return factory->CreateRawPtr(alias);													\
}																							\
																							\
template<>																					\
inline std::shared_ptr<BaseClassName>														\
CreateSharedPtr<BaseClassName>(const std::string &alias)									\
{																							\
	auto factory = GY_FACTORY_NAME(BaseClassName)::GetInstancePtr();						\
	return factory->CreateSharedPtr(alias);													\
}																							\
																							\
template<>																					\
inline std::unique_ptr<BaseClassName>														\
CreateInstanceUniquePtr<BaseClassName>(const std::string& alias)							\
{																							\
	auto factory = GY_FACTORY_NAME(BaseClassName)::GetInstancePtr();						\
	return factory->CreateUniquePtr(alias);													\
}																							\
																							\
template<>																					\
inline BaseClassName*																		\
CreatePlacementPtr<BaseClassName>(const std::string& alias, void* placement)				\
{																							\
	auto factory = GY_FACTORY_NAME(BaseClassName)::GetInstancePtr();						\
	return factory->CreatePlacementPtr(alias, placement);									\
}																							\
																							\


#define GY_IMPLEMENTATION_DEF(BaseClassName, ClassName, Alias)								\
class GY_IMPLEMENTATION_NAME(ClassName)														\
{																							\
	using CreateMethod = std::function<BaseClassName*()>;									\
public:																						\
	GY_IMPLEMENTATION_NAME(ClassName)()														\
	{																						\
		auto factory = GY_FACTORY_NAME(BaseClassName)::GetInstancePtr();					\
		factory->Register<ClassName>(Alias);												\
	}																						\
} Impementation_##BaseClassName##_##ClassName##_Instance;

GY_NAMESPACE_END